//! Stereological microstructure quantification.
//!
//! Implements the ASTM E112 mean-linear-intercept grain-size
//! method, area/volume fraction estimators from 2-D and 3-D image
//! data, and related point/line counting helpers.
//!
//! All functions operate on simple `f64` data so they can ingest
//! results from any image-analysis pipeline or phase-field solver.

use serde::{Deserialize, Serialize};

/// ASTM E112 mean-linear-intercept grain size.  Given a list of
/// horizontal line segments that intersect grain boundaries,
/// returns the mean intercept length `ℓ̄` and the ASTM grain-size
/// number `G`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LinearIntercept {
    /// Total line length scanned (m).
    pub total_line_length: f64,
    /// Number of grain-boundary intersections counted.
    pub n_intersections: usize,
}

impl LinearIntercept {
    /// Construct from a total line length and intersection count.
    pub fn new(total_line_length: f64, n_intersections: usize) -> Self {
        Self {
            total_line_length,
            n_intersections,
        }
    }

    /// Mean linear intercept `ℓ̄ = L / N` (m).
    pub fn mean_intercept(&self) -> f64 {
        if self.n_intersections == 0 {
            0.0
        } else {
            self.total_line_length / self.n_intersections as f64
        }
    }

    /// ASTM E112 grain-size number
    /// `G = -6.644 · log10(ℓ̄ / 1_inch_in_m)` (the negative sign is
    /// dropped to follow ASTM convention).
    pub fn astm_grain_size_number(&self) -> f64 {
        let l_in = self.mean_intercept() / 0.0254;
        if l_in <= 0.0 {
            0.0
        } else {
            -6.6439 * l_in.log10()
        }
    }
}

/// 2-D phase fraction from a binary pixel grid.  `n_phase_pixels`
/// is the number of pixels classified as belonging to the phase;
/// `n_total` is the total number of pixels in the field of view.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AreaFraction2D {
    /// Number of pixels in the phase.
    pub n_phase_pixels: usize,
    /// Total pixels.
    pub n_total: usize,
}

impl AreaFraction2D {
    /// Construct.
    pub fn new(n_phase_pixels: usize, n_total: usize) -> Self {
        Self {
            n_phase_pixels,
            n_total,
        }
    }

    /// `A_A = n_phase / n_total`.
    pub fn area_fraction(&self) -> f64 {
        if self.n_total == 0 {
            0.0
        } else {
            self.n_phase_pixels as f64 / self.n_total as f64
        }
    }

    /// Uncertainty from binomial counting `σ = sqrt(A_A (1 − A_A) / N)`.
    pub fn uncertainty(&self) -> f64 {
        let n = self.n_total as f64;
        if n == 0.0 {
            0.0
        } else {
            let a = self.area_fraction();
            (a * (1.0 - a) / n).sqrt()
        }
    }
}

/// 3-D volume fraction from voxel data.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VolumeFraction3D {
    /// Number of voxels in the phase.
    pub n_phase_voxels: usize,
    /// Total voxels.
    pub n_total: usize,
}

impl VolumeFraction3D {
    /// Construct.
    pub fn new(n_phase_voxels: usize, n_total: usize) -> Self {
        Self {
            n_phase_voxels,
            n_total,
        }
    }

    /// `V_V = n_phase / n_total`.
    pub fn volume_fraction(&self) -> f64 {
        if self.n_total == 0 {
            0.0
        } else {
            self.n_phase_voxels as f64 / self.n_total as f64
        }
    }
}

/// Mean number of features per unit area from a Delesse / point-
/// counting analysis.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NumberPerArea {
    /// Number of features counted.
    pub n_features: usize,
    /// Total scanned area (m²).
    pub total_area: f64,
}

impl NumberPerArea {
    /// `N_A = n_features / total_area` (m⁻²).
    pub fn number_per_area(&self) -> f64 {
        if self.total_area <= 0.0 {
            0.0
        } else {
            self.n_features as f64 / self.total_area
        }
    }
}

/// Estimate the 3-D number density from a 2-D section by the
/// Saltykov stereological method: given a histogram of section
/// circle radii with widths `Δr`, return the 3-D size-class
/// populations.
pub fn saltykov_size_distribution(
    radii_2d: &[f64],
    bin_width: f64,
) -> Vec<f64> {
    // Saltykov (1967): for each size class k,
    //   N_V(k) = (1 / ΔV_k) · Σ_{j=k}^{n_max} (-1)^{j-k} c(j,k) · N_A(j)
    // with `c(j,k)` the standard Saltykov coefficients and `ΔV_k`
    // the class-volume element.  We implement the simple uniform
    // case where the class volume is `4π r_k² Δr`.
    let mut counts: Vec<usize> = vec![0; radii_2d.len()];
    for &r in radii_2d {
        let idx = ((r / bin_width).floor() as usize).min(counts.len() - 1);
        counts[idx] += 1;
    }
    let mut result = vec![0.0_f64; counts.len()];
    for (k, &n_k) in counts.iter().enumerate() {
        let r_k = (k as f64 + 0.5) * bin_width;
        let d_v = 4.0 * std::f64::consts::PI * r_k.powi(2) * bin_width;
        if d_v > 0.0 {
            result[k] = n_k as f64 / d_v;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn linear_intercept_recovers_length() {
        let l = LinearIntercept::new(100.0e-6, 50);
        assert!(approx(l.mean_intercept(), 2.0e-6, 1.0e-12));
    }

    #[test]
    fn astm_grain_size_zero_intersections_returns_zero() {
        let l = LinearIntercept::new(1.0e-3, 0);
        assert_eq!(l.astm_grain_size_number(), 0.0);
    }

    #[test]
    fn area_fraction_matches_ratio() {
        let a = AreaFraction2D::new(30, 100);
        assert!(approx(a.area_fraction(), 0.30, 1.0e-9));
    }

    #[test]
    fn area_fraction_uncertainty_decreases_with_sample_size() {
        let a1 = AreaFraction2D::new(50, 100);
        let a2 = AreaFraction2D::new(500, 1000);
        assert!(a2.uncertainty() < a1.uncertainty());
    }

    #[test]
    fn volume_fraction_matches_ratio() {
        let v = VolumeFraction3D::new(700, 1000);
        assert!(approx(v.volume_fraction(), 0.7, 1.0e-9));
    }

    #[test]
    fn saltykov_size_distribution_non_negative() {
        let radii: Vec<f64> = vec![0.5, 1.5, 0.5, 2.5, 1.5];
        let dist = saltykov_size_distribution(&radii, 1.0);
        for v in &dist {
            assert!(*v >= 0.0);
        }
    }
}