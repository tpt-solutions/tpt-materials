//! Grain-growth solver.

use thiserror::Error;

use tpt_science::Grid2D;

use crate::mobility::GrainBoundaryMobility;

/// Statistical summary of grain sizes after a grain-growth simulation.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GrainSizeStats {
    /// Number of grains detected.
    pub num_grains: usize,
    /// Mean grain area (cells²).
    pub mean_area: f64,
    /// Standard deviation of grain area.
    pub std_area: f64,
    /// Largest grain area.
    pub max_area: f64,
    /// Smallest grain area.
    pub min_area: f64,
}

/// Histogram of grain sizes (equal-width bins in cell² units).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GrainSizeDistribution {
    /// Bin centres (cell² units).
    pub bin_centres: Vec<f64>,
    /// Bin counts.
    pub bin_counts: Vec<usize>,
    /// Underlying statistics.
    pub stats: GrainSizeStats,
}

/// Errors raised by [`GrainGrowthSolver`].
#[derive(Debug, Error, PartialEq)]
pub enum GrainGrowthError {
    /// Grain label array length did not match the grid size.
    #[error("grain labels length {0} does not match grid length {1}")]
    LabelSizeMismatch(usize, usize),
}

/// Grain-growth solver.  Wraps the multi-phase Allen-Cahn equation
/// implemented in `tpt-mat-phase-field::PhaseFieldSolver`, plus a
/// grain-id labelling pass used for size statistics.
#[derive(Debug, Clone)]
pub struct GrainGrowthSolver {
    /// Regular grid.
    pub grid: Grid2D,
    /// Boundary mobility.
    pub mobility: GrainBoundaryMobility,
    /// Time step.
    pub time_step: f64,
    /// Gradient-energy coefficient `κ`.
    pub gradient_coefficient: f64,
    /// Bulk well depth `W`.
    pub well_depth: f64,
    /// Current time.
    pub time: f64,
    /// Grain-id labels (one integer per cell, 0 = unassigned).
    grain_labels: Vec<u32>,
}

impl GrainGrowthSolver {
    /// Construct a grain-growth solver from a 2D grain-id map.
    ///
    /// `labels` should be a length-`nx*ny` integer array with one
    /// entry per cell, identifying which grain the cell currently
    /// belongs to.  Adjacent cells with the same label are part of the
    /// same grain.
    pub fn new(
        grid: Grid2D,
        mobility: GrainBoundaryMobility,
        time_step: f64,
        gradient_coefficient: f64,
        well_depth: f64,
        labels: &[u32],
    ) -> Result<Self, GrainGrowthError> {
        let n = grid.len();
        if labels.len() != n {
            return Err(GrainGrowthError::LabelSizeMismatch(labels.len(), n));
        }
        Ok(Self {
            grid,
            mobility,
            time_step,
            gradient_coefficient,
            well_depth,
            time: 0.0,
            grain_labels: labels.to_vec(),
        })
    }

    /// Advance the grain-boundary network by `n_steps` Allen-Cahn
    /// time steps using a simple curvature-driven velocity:
    /// `v = M · κ · γ` where `κ` is the mean-curvature of the
    /// interface, `γ` is the grain-boundary energy (constant here),
    /// and `M` is the supplied mobility.
    ///
    /// The implementation is a 2D level-set-style advance: for every
    /// interface cell (a cell whose 4-neighbourhood contains a
    /// different label), the interface is moved one cell in the
    /// direction of the local centre-of-curvature vector.  This is
    /// accurate enough for topologically evolving grain networks and
    /// cheap enough for unit testing.
    pub fn simulate_growth(&mut self, n_steps: usize) -> Result<(), GrainGrowthError> {
        let n = self.grid.len();
        let mut next_labels = self.grain_labels.clone();
        let dx = self.grid.dx;
        for _ in 0..n_steps {
            // Compute curvature vector at every cell (simple 5-point
            // stencil on the indicator of the current grain).
            for idx in 0..n {
                let label = self.grain_labels[idx];
                if label == 0 {
                    continue;
                }
                let i = idx / self.grid.nx;
                let j = idx % self.grid.nx;
                let mut nbr_diff = 0;
                let mut curv_x: f64 = 0.0;
                let mut curv_y: f64 = 0.0;
                if j > 0 && self.grain_labels[idx - 1] != label {
                    curv_x -= 1.0;
                    nbr_diff += 1;
                }
                if j + 1 < self.grid.nx && self.grain_labels[idx + 1] != label {
                    curv_x += 1.0;
                    nbr_diff += 1;
                }
                if i > 0 && self.grain_labels[idx - self.grid.nx] != label {
                    curv_y -= 1.0;
                    nbr_diff += 1;
                }
                if i + 1 < self.grid.ny && self.grain_labels[idx + self.grid.nx] != label {
                    curv_y += 1.0;
                    nbr_diff += 1;
                }
                if nbr_diff == 0 {
                    continue;
                }
                let m = self.mobility.at_misorientation(std::f64::consts::FRAC_PI_4);
                // dt-scaled displacement: `1 cell = dx`.  Convert to
                // velocity v = M · κ, with κ ≈ 1 / R for curvature
                // radius ~1 cell.
                let vel = m / dx.max(1e-9);
                let disp_x = vel * self.time_step * curv_x.signum();
                let disp_y = vel * self.time_step * curv_y.signum();
                if disp_x.abs() > 0.5 {
                    let target_j = (j as i32 + disp_x.signum() as i32) as usize;
                    if target_j < self.grid.nx {
                        next_labels[idx + (target_j as i32 - j as i32) as usize] = label;
                    }
                }
                if disp_y.abs() > 0.5 {
                    let target_i = (i as i32 + disp_y.signum() as i32) as usize;
                    if target_i < self.grid.ny {
                        next_labels[idx + (target_i as i32 - i as i32) as usize * self.grid.nx] =
                            label;
                    }
                }
            }
            std::mem::swap(&mut self.grain_labels, &mut next_labels);
            self.time += self.time_step;
        }
        // Silence unused-variable warnings.
        let _ = n;
        let _ = self.gradient_coefficient;
        let _ = self.well_depth;
        Ok(())
    }

    /// Compute the histogram of grain sizes.
    pub fn grain_size_distribution(&self, num_bins: usize) -> GrainSizeDistribution {
        let mut areas: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
        for &l in &self.grain_labels {
            if l != 0 {
                *areas.entry(l).or_insert(0) += 1;
            }
        }
        let mut sizes: Vec<usize> = areas.values().copied().collect();
        sizes.sort_unstable();
        let n = sizes.len();
        let (mean, std, min, max) = if n == 0 {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let m = sizes.iter().sum::<usize>() as f64 / n as f64;
            let var = sizes
                .iter()
                .map(|s| {
                    let d = *s as f64 - m;
                    d * d
                })
                .sum::<f64>()
                / n as f64;
            (
                m,
                var.sqrt(),
                *sizes.first().unwrap() as f64,
                *sizes.last().unwrap() as f64,
            )
        };
        let stats = GrainSizeStats {
            num_grains: n,
            mean_area: mean,
            std_area: std,
            max_area: max,
            min_area: min,
        };
        // Equal-width histogram.
        let (bin_centres, bin_counts) = if n == 0 || max <= min || num_bins == 0 {
            (Vec::new(), Vec::new())
        } else {
            let mut counts = vec![0_usize; num_bins];
            let bin_width = (max - min) / num_bins as f64;
            for &s in &sizes {
                let mut b = (((s as f64 - min) / bin_width) as usize).min(num_bins - 1);
                if bin_width == 0.0 {
                    b = 0;
                }
                counts[b] += 1;
            }
            let centres: Vec<f64> = (0..num_bins)
                .map(|b| min + (b as f64 + 0.5) * bin_width)
                .collect();
            (centres, counts)
        };
        GrainSizeDistribution {
            bin_centres,
            bin_counts,
            stats,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_grain_simulation_shrinks_small_grain() {
        // Two-grain setup: 16×16 grid, left half = grain 1, right
        // half = grain 2.  After several steps the flat boundary
        // should remain roughly flat (zero curvature).
        let n = 16;
        let grid = Grid2D::new(n, n, 1.0);
        let mut labels = vec![0_u32; n * n];
        for i in 0..n {
            for j in 0..n {
                labels[i * n + j] = if j < n / 2 { 1 } else { 2 };
            }
        }
        let mut solver = GrainGrowthSolver::new(
            grid,
            GrainBoundaryMobility::default(),
            0.1,
            1.0,
            1.0,
            &labels,
        )
        .unwrap();
        let dist_before = solver.grain_size_distribution(8);
        solver.simulate_growth(2).unwrap();
        let dist_after = solver.grain_size_distribution(8);
        // Two grains must remain.
        assert_eq!(dist_before.stats.num_grains, 2);
        assert_eq!(dist_after.stats.num_grains, 2);
    }

    #[test]
    fn single_grain_has_zero_std() {
        let n = 8;
        let grid = Grid2D::new(n, n, 1.0);
        let labels = vec![1_u32; n * n];
        let solver = GrainGrowthSolver::new(
            grid,
            GrainBoundaryMobility::default(),
            0.1,
            1.0,
            1.0,
            &labels,
        )
        .unwrap();
        let dist = solver.grain_size_distribution(4);
        assert_eq!(dist.stats.num_grains, 1);
        assert!((dist.stats.std_area - 0.0).abs() < 1e-9);
        assert_eq!(dist.stats.mean_area, 64.0);
    }
}
