//! Phase-field solver: time integration of Allen-Cahn, Cahn-Hilliard,
//! Kobayashi, and multi-phase models on regular grids.

use thiserror::Error;

use tpt_science::{Grid2D, Grid3D};

use crate::energy::{BulkEnergy, FreeEnergyFunctional};
use crate::model::PhaseFieldModel;
use crate::result::PhaseFieldResult;

/// Errors raised by [`PhaseFieldSolver`].
#[derive(Debug, Error, PartialEq)]
pub enum PhaseFieldError {
    /// Field size did not match the grid size.
    #[error("field length {0} does not match grid length {1}")]
    FieldSizeMismatch(usize, usize),
    /// Unknown model configuration.
    #[error("unknown model configuration")]
    UnknownModel,
}

/// Phase-field solver configuration.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseFieldSolver {
    /// Model (Allen-Cahn, Cahn-Hilliard, Kobayashi, MultiPhase).
    pub model: PhaseFieldModel,
    /// 2D grid (square cells).
    pub grid_2d: Grid2D,
    /// Time step.
    pub time_step: f64,
    /// Mobility / kinetic coefficient `L` (Allen-Cahn, multi-phase).
    pub mobility: f64,
    /// Gradient-energy coefficient `κ`.
    pub gradient_coefficient: f64,
    /// Diffusion coefficient (Cahn-Hilliard mobility `M`).
    pub diffusivity: f64,
    /// Latent-heat coupling (Kobayashi).
    pub latent_heat: f64,
    /// Thermal diffusivity (Kobayashi).
    pub thermal_diffusivity: f64,
    /// Bulk free energy.
    pub bulk_energy: BulkEnergy,
    /// Current simulation time.
    pub time: f64,
    /// Order-parameter field(s).
    fields: Vec<Vec<f64>>,
    /// Concentration (Cahn-Hilliard only).
    pub concentration: Vec<f64>,
    /// Temperature (Kobayashi only).
    pub temperature: Vec<f64>,
}

impl PhaseFieldSolver {
    /// Set the temperature field (Kobayashi solver).
    pub fn set_temperature(&mut self, t: Vec<f64>) {
        if t.len() == self.temperature.len() {
            self.temperature = t;
        }
    }

    /// Mutable access to the temperature field (Kobayashi solver).
    pub fn temperature_mut(&mut self) -> &mut [f64] {
        &mut self.temperature
    }

    /// Build a solver for a Cahn-Hilliard simulation, where the input is
    /// an initial concentration field.
    pub fn new_cahn_hilliard(
        grid: Grid2D,
        time_step: f64,
        diffusivity: f64,
        gradient_coefficient: f64,
        bulk_energy: BulkEnergy,
        initial_concentration: &[f64],
    ) -> Result<Self, PhaseFieldError> {
        let n = grid.len();
        if initial_concentration.len() != n {
            return Err(PhaseFieldError::FieldSizeMismatch(
                initial_concentration.len(),
                n,
            ));
        }
        Ok(Self {
            model: PhaseFieldModel::CahnHilliard,
            grid_2d: grid,
            time_step,
            mobility: 1.0,
            gradient_coefficient,
            diffusivity,
            latent_heat: 0.0,
            thermal_diffusivity: 0.0,
            bulk_energy,
            time: 0.0,
            fields: vec![initial_concentration.to_vec()],
            concentration: initial_concentration.to_vec(),
            temperature: vec![0.0; n],
        })
    }
    pub fn new_2d(
        model: PhaseFieldModel,
        grid: Grid2D,
        time_step: f64,
        mobility: f64,
        gradient_coefficient: f64,
        bulk_energy: BulkEnergy,
        initial_order: &[f64],
    ) -> Result<Self, PhaseFieldError> {
        let n = grid.len();
        if initial_order.len() != n {
            return Err(PhaseFieldError::FieldSizeMismatch(initial_order.len(), n));
        }
        Ok(Self {
            model,
            grid_2d: grid,
            time_step,
            mobility,
            gradient_coefficient,
            diffusivity: 0.0,
            latent_heat: 0.0,
            thermal_diffusivity: 0.0,
            bulk_energy,
            time: 0.0,
            fields: vec![initial_order.to_vec()],
            concentration: vec![0.5; n],
            temperature: vec![0.0; n],
        })
    }

    /// Number of grid cells.
    pub fn n_cells(&self) -> usize {
        self.grid_2d.len()
    }

    /// Current order parameter field (length 1 for single-phase
    /// models, `num_grains` for multi-phase).
    pub fn order_parameter(&self) -> &[Vec<f64>] {
        &self.fields
    }

    /// Advance the simulation by `n_steps` time steps.
    pub fn step_many(&mut self, n_steps: usize) -> Result<(), PhaseFieldError> {
        for _ in 0..n_steps {
            self.step()?;
        }
        Ok(())
    }

    /// Advance by exactly one time step.
    pub fn step(&mut self) -> Result<(), PhaseFieldError> {
        match self.model {
            PhaseFieldModel::AllenCahn => self.step_allen_cahn(),
            PhaseFieldModel::CahnHilliard => self.step_cahn_hilliard(),
            PhaseFieldModel::Kobayashi => self.step_kobayashi(),
            PhaseFieldModel::MultiPhase { num_grains } => self.step_multi_phase(num_grains),
        }
    }

    /// One Allen-Cahn step:
    /// `η^{n+1} = η^n + Δt · L (κ ∇²η − f′(η))`.
    pub fn step_allen_cahn(&mut self) -> Result<(), PhaseFieldError> {
        let grid = self.grid_2d.clone();
        let n = grid.len();
        let field = &self.fields[0];
        if field.len() != n {
            return Err(PhaseFieldError::FieldSizeMismatch(field.len(), n));
        }
        let mut lap = vec![0.0_f64; n];
        grid.laplacian_neumann(field, &mut lap);
        let mut next = field.to_vec();
        let dt = self.time_step;
        let kappa = self.gradient_coefficient;
        let l = self.mobility;
        let bulk = self.bulk_energy;
        for i in 0..n {
            let df = bulk.derivative(field[i]);
            next[i] = field[i] + dt * l * (kappa * lap[i] - df);
        }
        self.fields[0] = next;
        self.time += dt;
        Ok(())
    }

    /// One Cahn-Hilliard step.
    pub fn step_cahn_hilliard(&mut self) -> Result<(), PhaseFieldError> {
        let grid = self.grid_2d.clone();
        let n = grid.len();
        let c = self.concentration.clone();
        if c.len() != n {
            return Err(PhaseFieldError::FieldSizeMismatch(c.len(), n));
        }
        let mut lap_c = vec![0.0_f64; n];
        grid.laplacian_neumann(&c, &mut lap_c);
        let bulk = self.bulk_energy;
        let mut mu = vec![0.0_f64; n];
        for i in 0..n {
            mu[i] = bulk.derivative(c[i]) - self.gradient_coefficient * lap_c[i];
        }
        let mut j_x = vec![0.0_f64; n];
        let mut j_y = vec![0.0_f64; n];
        let dx = grid.dx;
        for i in 0..grid.ny {
            for j in 0..grid.nx {
                let c_idx = grid.idx(i, j);
                let mu_right = if j + 1 < grid.nx {
                    mu[grid.idx(i, j + 1)]
                } else {
                    mu[grid.idx(i, j - 1)]
                };
                let mu_left = if j > 0 {
                    mu[grid.idx(i, j - 1)]
                } else {
                    mu[grid.idx(i, j + 1)]
                };
                let mu_up = if i + 1 < grid.ny {
                    mu[grid.idx(i + 1, j)]
                } else {
                    mu[grid.idx(i - 1, j)]
                };
                let mu_down = if i > 0 {
                    mu[grid.idx(i - 1, j)]
                } else {
                    mu[grid.idx(i + 1, j)]
                };
                j_x[c_idx] = -self.diffusivity * (mu_right - mu_left) * 0.5 / dx;
                j_y[c_idx] = -self.diffusivity * (mu_up - mu_down) * 0.5 / dx;
            }
        }
        let mut next = vec![0.0_f64; n];
        for i in 0..grid.ny {
            for j in 0..grid.nx {
                let idx = grid.idx(i, j);
                let j_right = if j + 1 < grid.nx {
                    j_x[grid.idx(i, j + 1)]
                } else {
                    j_x[grid.idx(i, j - 1)]
                };
                let j_left = if j > 0 {
                    j_x[grid.idx(i, j - 1)]
                } else {
                    j_x[grid.idx(i, j + 1)]
                };
                let j_up = if i + 1 < grid.ny {
                    j_y[grid.idx(i + 1, j)]
                } else {
                    j_y[grid.idx(i - 1, j)]
                };
                let j_down = if i > 0 {
                    j_y[grid.idx(i - 1, j)]
                } else {
                    j_y[grid.idx(i + 1, j)]
                };
                let div = (j_right - j_left) * 0.5 / dx + (j_up - j_down) * 0.5 / dx;
                next[idx] = c[idx] + self.time_step * (-div);
            }
        }
        self.concentration = next;
        if self.fields.is_empty() {
            self.fields.push(self.concentration.clone());
        } else {
            self.fields[0] = self.concentration.clone();
        }
        self.time += self.time_step;
        Ok(())
    }

    /// One Kobayashi step.
    pub fn step_kobayashi(&mut self) -> Result<(), PhaseFieldError> {
        let grid = self.grid_2d.clone();
        let n = grid.len();
        let eta = self.fields[0].clone();
        let u = self.temperature.clone();
        let mut lap_eta = vec![0.0_f64; n];
        let mut lap_u = vec![0.0_f64; n];
        grid.laplacian_neumann(&eta, &mut lap_eta);
        grid.laplacian_neumann(&u, &mut lap_u);
        let dt = self.time_step;
        let bulk = self.bulk_energy;
        let mut new_eta = eta.clone();
        for i in 0..n {
            let df = bulk.derivative(eta[i]);
            new_eta[i] = eta[i]
                + dt * self.mobility
                    * (self.gradient_coefficient * lap_eta[i] - df + self.latent_heat * u[i]);
        }
        let mut new_u = vec![0.0_f64; n];
        for i in 0..n {
            new_u[i] = u[i]
                + dt * (self.thermal_diffusivity * lap_u[i]
                    + 0.5 * self.mobility * (new_eta[i] - eta[i]));
        }
        self.fields[0] = new_eta;
        self.temperature = new_u;
        self.time += dt;
        Ok(())
    }

    /// One multi-phase Allen-Cahn step (simplified Moelans version).
    pub fn step_multi_phase(&mut self, num_grains: usize) -> Result<(), PhaseFieldError> {
        if self.fields.len() != num_grains {
            return Err(PhaseFieldError::UnknownModel);
        }
        let grid = self.grid_2d.clone();
        let n = grid.len();
        let dt = self.time_step;
        let kappa = self.gradient_coefficient;
        let l = self.mobility;
        let mut next = self.fields.clone();
        let mut lap_buf = vec![0.0_f64; n];
        for (g, eta) in self.fields.iter().enumerate() {
            grid.laplacian_neumann(eta, &mut lap_buf);
            for i in 0..n {
                let df = self.bulk_energy.derivative(eta[i]);
                let mut coupling = 0.0;
                for (h, other) in self.fields.iter().enumerate() {
                    if h != g {
                        coupling += other[i] * other[i];
                    }
                }
                next[g][i] = eta[i] + dt * l * (kappa * lap_buf[i] - df + 2.0 * eta[i] * coupling);
            }
        }
        self.fields = next;
        self.time += dt;
        Ok(())
    }

    /// Build a [`PhaseFieldResult`] snapshot at the current time.
    pub fn snapshot(&self) -> PhaseFieldResult {
        let free_energy = if matches!(self.model, PhaseFieldModel::CahnHilliard) {
            let fe = FreeEnergyFunctional {
                gradient_coefficient: self.gradient_coefficient,
                dx: self.grid_2d.dx,
            };
            fe.evaluate_2d(
                &self.concentration,
                self.grid_2d.nx,
                self.grid_2d.ny,
                &self.bulk_energy,
            )
        } else {
            let fe = FreeEnergyFunctional {
                gradient_coefficient: self.gradient_coefficient,
                dx: self.grid_2d.dx,
            };
            fe.evaluate_2d(
                &self.fields[0],
                self.grid_2d.nx,
                self.grid_2d.ny,
                &self.bulk_energy,
            )
        };
        let threshold = 0.05;
        let mut area = 0.0;
        let field = &self.fields[0];
        let nx = self.grid_2d.nx;
        let ny = self.grid_2d.ny;
        let dx = self.grid_2d.dx;
        for i in 0..ny {
            for j in 0..nx {
                let c = i * nx + j;
                let er = if j + 1 < nx {
                    field[c + 1]
                } else {
                    field[c - 1]
                };
                let el = if j > 0 { field[c - 1] } else { field[c + 1] };
                let eu = if i + 1 < ny {
                    field[c + nx]
                } else {
                    field[c - nx]
                };
                let ed = if i > 0 { field[c - nx] } else { field[c + nx] };
                let gx = (er - el) * 0.5 / dx;
                let gy = (eu - ed) * 0.5 / dx;
                let grad_norm = (gx * gx + gy * gy).sqrt();
                if grad_norm > threshold {
                    area += 1.0;
                }
            }
        }
        area *= dx * dx;
        PhaseFieldResult {
            order_parameter: self.fields.clone(),
            concentration: self.concentration.clone(),
            temperature: self.temperature.clone(),
            free_energy,
            interface_area: area,
            time: self.time,
        }
    }

    /// Stand-alone helper: compute the discrete Laplacian of a 2D
    /// field via Neumann zero-flux boundaries.
    pub fn compute_laplacian(field: &[f64], grid: &Grid2D, out: &mut [f64]) {
        grid.laplacian_neumann(field, out);
    }

    /// Apply a per-cell multiplicative correction to the first order
    /// parameter field.  Used by `tpt-mat-solidification` to encode
    /// anisotropy on the interface cells.
    pub fn apply_field_correction(&mut self, f: impl Fn(f64) -> f64) {
        for v in self.fields[0].iter_mut() {
            *v = f(*v);
        }
    }

    /// Read-only access to the grid.
    pub fn grid(&self) -> &Grid2D {
        &self.grid_2d
    }

    /// Convenience for `Grid3D` users: explicit 3D Laplacian on the
    /// supplied field.  Lives here so the phase-field API stays the
    /// single entry point for `tpt-science::grid` consumers.
    pub fn compute_laplacian_3d(field: &[f64], grid: &Grid3D, out: &mut [f64]) {
        grid.laplacian_neumann(field, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn small_grid(n: usize) -> Grid2D {
        Grid2D::new(n, n, 1.0)
    }

    #[test]
    fn allen_cahn_decreases_free_energy() {
        let n = 32;
        let grid = small_grid(n);
        let mut eta = vec![0.5_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                eta[i * n + j] = if j < n / 2 { 0.0 } else { 1.0 };
            }
        }
        let mut solver = PhaseFieldSolver::new_2d(
            PhaseFieldModel::AllenCahn,
            grid,
            0.05,
            1.0,
            1.0,
            BulkEnergy::DoubleWell { well_depth: 1.0 },
            &eta,
        )
        .unwrap();
        let e0 = solver.snapshot().free_energy;
        for _ in 0..50 {
            solver.step().unwrap();
        }
        let e1 = solver.snapshot().free_energy;
        assert!(e1 <= e0 + 1e-9, "energy increased: {e0} → {e1}");
    }

    #[test]
    fn cahn_hilliard_conserves_average_concentration() {
        let n = 16;
        let grid = small_grid(n);
        let mut c = vec![0.5; n * n];
        for (i, v) in c.iter_mut().enumerate() {
            let x = (i % n) as f64 / n as f64;
            *v = 0.5 + 0.05 * (2.0 * std::f64::consts::PI * x).sin();
        }
        let mean0: f64 = c.iter().sum::<f64>() / c.len() as f64;
        let mut solver = PhaseFieldSolver::new_2d(
            PhaseFieldModel::CahnHilliard,
            grid,
            0.001,
            1.0,
            1.0,
            BulkEnergy::RegularSolution(crate::energy::RegularSolutionParams {
                omega: 5.0,
                rt: 1.0,
            }),
            &c,
        )
        .unwrap();
        solver.diffusivity = 1.0;
        for _ in 0..20 {
            solver.step().unwrap();
        }
        let mean1: f64 =
            solver.concentration.iter().sum::<f64>() / solver.concentration.len() as f64;
        assert_relative_eq!(mean0, mean1, epsilon = 1e-9);
    }

    #[test]
    fn kobayashi_zero_latent_heat_matches_allen_cahn() {
        let n = 16;
        let grid = small_grid(n);
        let mut eta = vec![0.5; n * n];
        for (i, e) in eta.iter_mut().enumerate() {
            let x = i % n;
            *e = if x < n / 2 { 0.0 } else { 1.0 };
        }
        let mut solver = PhaseFieldSolver::new_2d(
            PhaseFieldModel::Kobayashi,
            grid,
            0.05,
            1.0,
            1.0,
            BulkEnergy::DoubleWell { well_depth: 1.0 },
            &eta,
        )
        .unwrap();
        solver.latent_heat = 0.0;
        solver.thermal_diffusivity = 1.0;
        for _ in 0..10 {
            solver.step().unwrap();
        }
        let g = solver.fields[0].clone();
        assert!(g.iter().any(|v| *v > 0.0 && *v < 1.0));
    }

    #[test]
    fn compute_laplacian_matches_grid_helper() {
        let grid = small_grid(8);
        let f: Vec<f64> = (0..64).map(|i| (i as f64) * 0.01).collect();
        let mut out = vec![0.0_f64; 64];
        PhaseFieldSolver::compute_laplacian(&f, &grid, &mut out);
        let mut out_ref = vec![0.0_f64; 64];
        grid.laplacian_neumann(&f, &mut out_ref);
        for (a, b) in out.iter().zip(out_ref.iter()) {
            assert_eq!(*a, *b);
        }
    }

    #[test]
    fn temperature_mut_seeds_undercooling() {
        let n = 8;
        let grid = small_grid(n);
        let eta = vec![0.0_f64; n * n];
        let mut solver = PhaseFieldSolver::new_2d(
            PhaseFieldModel::Kobayashi,
            grid,
            0.01,
            1.0,
            1.0,
            BulkEnergy::DoubleWell { well_depth: 1.0 },
            &eta,
        )
        .unwrap();
        for u in solver.temperature_mut() {
            *u = -1.0;
        }
        assert!(solver.temperature.iter().all(|&u| u == -1.0));
    }
}
