//! Single-component Fick's-second-law solver on regular grids.

use thiserror::Error;

use tpt_science::Grid2D;

/// Errors raised by [`DiffusionSolver`].
#[derive(Debug, Error, PartialEq)]
pub enum DiffusionSolverError {
    /// Concentration field size did not match the grid size.
    #[error("concentration length {0} does not match grid length {1}")]
    ConcentrationSizeMismatch(usize, usize),
    /// Negative time step supplied.
    #[error("time step must be positive, got {0}")]
    NegativeTimeStep(f64),
    /// Negative diffusivity supplied.
    #[error("diffusivity must be non-negative, got {0}")]
    NegativeDiffusivity(f64),
}

/// Result of a single diffusion step.
#[derive(Debug, Clone, PartialEq)]
pub struct DiffusionStepResult {
    /// Number of cells.
    pub n_cells: usize,
    /// Sum of concentration (should be conserved).
    pub total_mass: f64,
    /// Maximum absolute concentration change in the step.
    pub max_abs_change: f64,
}

/// Single-component diffusion solver.
///
/// `∂c/∂t = D ∇²c`
///
/// Spatial discretisation: 5-point central-difference Laplacian with
/// Neumann zero-flux boundaries (`tpt_science::Grid2D::laplacian_neumann`).
/// Temporal discretisation: forward Euler, with the textbook stability
/// limit `Δt ≤ Δx² / (4 D)`.
#[derive(Debug, Clone)]
pub struct DiffusionSolver {
    /// Grid.
    pub grid: Grid2D,
    /// Diffusivity `D` (constant per call).
    pub diffusivity: f64,
    /// Time step.
    pub time_step: f64,
    /// Concentration field.
    pub concentration: Vec<f64>,
    /// Current time.
    pub time: f64,
}

impl DiffusionSolver {
    /// Construct a new solver from an initial concentration field.
    pub fn new(
        grid: Grid2D,
        diffusivity: f64,
        time_step: f64,
        initial_concentration: &[f64],
    ) -> Result<Self, DiffusionSolverError> {
        if diffusivity < 0.0 {
            return Err(DiffusionSolverError::NegativeDiffusivity(diffusivity));
        }
        if time_step <= 0.0 {
            return Err(DiffusionSolverError::NegativeTimeStep(time_step));
        }
        let n = grid.len();
        if initial_concentration.len() != n {
            return Err(DiffusionSolverError::ConcentrationSizeMismatch(
                initial_concentration.len(),
                n,
            ));
        }
        Ok(Self {
            grid,
            diffusivity,
            time_step,
            concentration: initial_concentration.to_vec(),
            time: 0.0,
        })
    }

    /// CFL stability limit `Δt ≤ Δx² / (4 D)` for forward Euler.
    pub fn cfl_limit(&self) -> f64 {
        let dx2 = self.grid.dx * self.grid.dx;
        if self.diffusivity <= 0.0 {
            f64::INFINITY
        } else {
            dx2 / (4.0 * self.diffusivity)
        }
    }

    /// Advance the concentration by one time step.
    pub fn step(&mut self) -> DiffusionStepResult {
        let n = self.grid.len();
        let mut lap = vec![0.0_f64; n];
        self.grid.laplacian_neumann(&self.concentration, &mut lap);
        let mut next = self.concentration.clone();
        let dt = self.time_step;
        let d = self.diffusivity;
        let mut max_change: f64 = 0.0;
        for i in 0..n {
            next[i] = self.concentration[i] + dt * d * lap[i];
            let dc = (next[i] - self.concentration[i]).abs();
            if dc > max_change {
                max_change = dc;
            }
        }
        let total: f64 = next.iter().sum();
        self.concentration = next;
        self.time += dt;
        DiffusionStepResult {
            n_cells: n,
            total_mass: total,
            max_abs_change: max_change,
        }
    }

    /// Advance by `n_steps` time steps.
    pub fn step_many(&mut self, n_steps: usize) -> DiffusionStepResult {
        let mut last = DiffusionStepResult {
            n_cells: self.grid.len(),
            total_mass: self.concentration.iter().sum(),
            max_abs_change: 0.0,
        };
        for _ in 0..n_steps {
            last = self.step();
        }
        last
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_diffuses_outward_and_conserves_mass() {
        // 2D Gaussian initial profile.  Neumann zero-flux boundary
        // → total mass should be exactly preserved.
        let n = 32;
        let grid = Grid2D::new(n, n, 1.0);
        let mut c = vec![0.0_f64; n * n];
        let centre = (n as f64 - 1.0) * 0.5;
        let sigma2 = 2.0_f64;
        let mut mass = 0.0;
        for i in 0..n {
            for j in 0..n {
                let dx = i as f64 - centre;
                let dy = j as f64 - centre;
                let v = (-(dx * dx + dy * dy) / (2.0 * sigma2)).exp();
                c[i * n + j] = v;
                mass += v;
            }
        }
        let mut solver = DiffusionSolver::new(grid, 0.2, 0.5, &c).unwrap();
        let dt_max = solver.cfl_limit();
        assert!(solver.time_step < dt_max);
        let r0 = solver.step();
        let r1 = solver.step_many(20);
        // Mass conservation to ~relative 1e-6 (Neumann).
        assert!((r0.total_mass - mass).abs() / mass < 1e-6);
        assert!((r1.total_mass - mass).abs() / mass < 1e-6);
        // Peak should decrease monotonically over many steps.
        let peak: f64 = solver
            .concentration
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(peak < 1.0, "peak = {peak} (expected < 1)");
    }

    #[test]
    fn zero_diffusivity_is_identity() {
        let n = 8;
        let grid = Grid2D::new(n, n, 1.0);
        let c: Vec<f64> = (0..n * n).map(|i| (i as f64) * 0.01).collect();
        let mut solver = DiffusionSolver::new(grid, 0.0, 0.1, &c).unwrap();
        let r = solver.step();
        for (a, b) in solver.concentration.iter().zip(c.iter()) {
            assert_eq!(*a, *b);
        }
        assert_eq!(r.max_abs_change, 0.0);
    }

    #[test]
    fn negative_diffusivity_rejected() {
        let grid = Grid2D::new(4, 4, 1.0);
        let c = vec![0.0_f64; 16];
        assert!(DiffusionSolver::new(grid, -1.0, 0.1, &c).is_err());
    }

    #[test]
    fn wrong_field_size_rejected() {
        let grid = Grid2D::new(4, 4, 1.0);
        let c = vec![0.0_f64; 15];
        assert!(DiffusionSolver::new(grid, 1.0, 0.1, &c).is_err());
    }
}
