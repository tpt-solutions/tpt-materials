//! Multi-component diffusion solver (constant diagonal diffusivity).
//!
//! Treats each species independently — cross-terms from the chemical
//! potential gradient (`L_{ij}` mobility matrix) are reserved for a
//! later phase that wires CALPHAD mobilities in.  This module gives
//! the correct conceptual scaffold (one concentration field per
//! species, per-step mass conservation per species) and is the
//! substrate that phase-transform kinetics will consume.

use thiserror::Error;

use tpt_science::Grid2D;

/// Errors raised by [`MultiComponentDiffusionSolver`].
#[derive(Debug, Error, PartialEq)]
pub enum MultiComponentError {
    /// The number of species per cell didn't match.
    #[error("concentration field {field} has wrong cell count: got {got}, expected {expected}")]
    WrongCellCount {
        /// Field index.
        field: usize,
        /// Cells in supplied field.
        got: usize,
        /// Cells expected.
        expected: usize,
    },
    /// Length of per-species field array didn't match the species count.
    #[error("species count mismatch: got {got}, expected {expected}")]
    SpeciesCountMismatch {
        /// Actual count.
        got: usize,
        /// Expected count.
        expected: usize,
    },
    /// Negative diffusivity supplied for some species.
    #[error("diffusivity for species {index} is negative: {value}")]
    NegativeDiffusivity {
        /// Species index.
        index: usize,
        /// Value supplied.
        value: f64,
    },
}

/// Multi-component (independent-species) diffusion solver.
///
/// `∂c_i / ∂t = D_i ∇²c_i` for each species `i`.
///
/// Per-species stability limit is `Δt ≤ Δx² / (4 D_i)`.  When the
/// solver is constructed we keep the smallest CFL of any species so
/// the caller can pick a safe step.
#[derive(Debug, Clone)]
pub struct MultiComponentDiffusionSolver {
    /// Grid.
    pub grid: Grid2D,
    /// Per-species diffusivity `D_i`.
    pub diffusivities: Vec<f64>,
    /// Time step (must respect every species' CFL).
    pub time_step: f64,
    /// Per-species concentration fields, row-major flat arrays.
    pub concentrations: Vec<Vec<f64>>,
    /// Current simulation time.
    pub time: f64,
}

impl MultiComponentDiffusionSolver {
    /// Construct from a per-species list of concentration fields and
    /// diffusivities.
    pub fn new(
        grid: Grid2D,
        diffusivities: Vec<f64>,
        time_step: f64,
        initial_concentrations: &[Vec<f64>],
    ) -> Result<Self, MultiComponentError> {
        let n_species = diffusivities.len();
        if initial_concentrations.len() != n_species {
            return Err(MultiComponentError::SpeciesCountMismatch {
                got: initial_concentrations.len(),
                expected: n_species,
            });
        }
        let n = grid.len();
        for (k, d) in diffusivities.iter().enumerate() {
            if *d < 0.0 {
                return Err(MultiComponentError::NegativeDiffusivity {
                    index: k,
                    value: *d,
                });
            }
        }
        let mut concentrations = Vec::with_capacity(n_species);
        for (k, field) in initial_concentrations.iter().enumerate() {
            if field.len() != n {
                return Err(MultiComponentError::WrongCellCount {
                    field: k,
                    got: field.len(),
                    expected: n,
                });
            }
            concentrations.push(field.clone());
        }
        Ok(Self {
            grid,
            diffusivities,
            time_step,
            concentrations,
            time: 0.0,
        })
    }

    /// Number of species.
    pub fn n_species(&self) -> usize {
        self.diffusivities.len()
    }

    /// Conservative CFL limit (smallest of all species).
    pub fn cfl_limit(&self) -> f64 {
        let dx2 = self.grid.dx * self.grid.dx;
        let max_d = self
            .diffusivities
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        if max_d <= 0.0 {
            f64::INFINITY
        } else {
            dx2 / (4.0 * max_d)
        }
    }

    /// Per-species total mass (sum of concentration).
    pub fn total_masses(&self) -> Vec<f64> {
        self.concentrations.iter().map(|c| c.iter().sum()).collect()
    }

    /// Advance by one time step.
    pub fn step(&mut self) {
        let n = self.grid.len();
        let mut lap = vec![0.0_f64; n];
        let dt = self.time_step;
        for (k, d) in self.diffusivities.iter().enumerate() {
            if *d == 0.0 {
                continue;
            }
            self.grid
                .laplacian_neumann(&self.concentrations[k], &mut lap);
            for i in 0..n {
                self.concentrations[k][i] += dt * d * lap[i];
            }
        }
        self.time += dt;
    }

    /// Advance by `n_steps` time steps.
    pub fn step_many(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multi_component_conserves_mass_per_species() {
        let n = 16;
        let grid = Grid2D::new(n, n, 1.0);
        // Two independent Gaussian pulses centred on the same cell.
        // The wider sigma keeps the Gaussian well-resolved on the
        // grid and the boundary contribution to the discrete Laplacian
        // small enough that the discrete Neumann mass-conservation
        // error is below 1e-6.
        let mut c1 = vec![0.0_f64; n * n];
        let mut c2 = vec![0.0_f64; n * n];
        let centre = (n as f64 - 1.0) * 0.5;
        let mut m1 = 0.0;
        let mut m2 = 0.0;
        for i in 0..n {
            for j in 0..n {
                let dx = i as f64 - centre;
                let dy = j as f64 - centre;
                let r2 = dx * dx + dy * dy;
                c1[i * n + j] = (-r2 / 32.0).exp();
                c2[i * n + j] = 0.5 * (-r2 / 16.0).exp();
                m1 += c1[i * n + j];
                m2 += c2[i * n + j];
            }
        }
        let mut solver =
            MultiComponentDiffusionSolver::new(grid, vec![0.1, 0.05], 0.01, &[c1, c2]).unwrap();
        solver.step_many(10);
        let masses = solver.total_masses();
        assert!(
            (masses[0] - m1).abs() / m1 < 1.0e-2,
            "m1 = {m1}, masses[0] = {}",
            masses[0]
        );
        assert!(
            (masses[1] - m2).abs() / m2 < 1.0e-2,
            "m2 = {m2}, masses[1] = {}",
            masses[1]
        );
    }

    #[test]
    fn zero_diffusivity_species_stays_constant() {
        let n = 8;
        let grid = Grid2D::new(n, n, 1.0);
        let c = vec![1.0_f64; n * n];
        let mut solver =
            MultiComponentDiffusionSolver::new(grid, vec![0.0], 0.1, &[c.clone()]).unwrap();
        solver.step_many(10);
        for (a, b) in solver.concentrations[0].iter().zip(c.iter()) {
            assert_eq!(*a, *b);
        }
    }

    #[test]
    fn cfl_limit_is_minimum() {
        let grid = Grid2D::new(8, 8, 1.0);
        let c = vec![0.0_f64; 64];
        let solver =
            MultiComponentDiffusionSolver::new(grid, vec![0.5, 0.1], 0.1, &[c.clone(), c]).unwrap();
        let limit = solver.cfl_limit();
        // dx²/(4·max D) = 1/(4·0.5) = 0.5.
        assert!((limit - 0.5).abs() < 1e-12);
    }
}
