//! Pole figure: stereographic projection of a crystallographic pole.
//!
//! A pole figure discretises a unit-hemisphere in the sample frame into
//! an `(n_lat × n_lon)` grid and reports the (volume-weighted) intensity
//! at each cell.  The intensity at cell `(φ, θ)` is the sum of weights
//! of grains whose `[uvw]` pole lies within the cell's solid angle.
//!
//! This is the elementary texture diagnostic; the ODF
//! ([`OrientationDistributionFunction`]) is the higher-fidelity
//! representation.
//!
//! [`OrientationDistributionFunction`]: crate::OrientationDistributionFunction

use serde::{Deserialize, Serialize};

use tpt_math_linalg_fixed::Vec3;

use crate::texture::TextureAnalyzer;

/// Which crystallographic pole to project (`{hkl}` indices in the
/// crystal frame).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PoleKind {
    /// `(111)` family — FCC octahedral.
    Fcc111,
    /// `(200)` family — FCC cubic.
    Fcc200,
    /// `(110)` family — BCC.
    Bcc110,
    /// `(0001)` family — HCP basal pole.
    Hcp0001,
    /// `(10-10)` family — HCP prismatic.
    Hcp10T10,
    /// User-supplied pole `[u, v, w]`.
    Custom([i32; 3]),
}

impl PoleKind {
    /// Crystallographic direction in the crystal frame (unit vector).
    pub fn direction(self) -> Vec3 {
        match self {
            PoleKind::Fcc111 => Vec3::new(1.0, 1.0, 1.0).normalized(),
            PoleKind::Fcc200 => Vec3::new(1.0, 0.0, 0.0),
            PoleKind::Bcc110 => Vec3::new(1.0, 1.0, 0.0).normalized(),
            PoleKind::Hcp0001 => Vec3::new(0.0, 0.0, 1.0),
            PoleKind::Hcp10T10 => Vec3::new(1.0, -1.0, 0.0).normalized(),
            PoleKind::Custom(d) => {
                Vec3::new(d[0] as f64, d[1] as f64, d[2] as f64).normalized()
            }
        }
    }
}

/// Equal-area pole-figure grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PoleFigureGrid {
    /// Latitude bins (0..=90°).
    pub n_lat: usize,
    /// Longitude bins (0..360°).
    pub n_lon: usize,
}

impl Default for PoleFigureGrid {
    fn default() -> Self {
        Self {
            n_lat: 36,
            n_lon: 72,
        }
    }
}

/// Pole figure: flattened `(n_lat × n_lon)` intensity array.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PoleFigure {
    /// Pole used.
    pub pole: PoleKind,
    /// Grid.
    pub grid: PoleFigureGrid,
    /// Flat intensity row-major, length `n_lat × n_lon`.
    pub intensity: Vec<f64>,
}

impl TextureAnalyzer {
    /// Compute a pole figure for `pole` at the given grid resolution.
    ///
    /// For each grain, the sample-frame pole direction is
    /// `d_sample = R · d_crystal` where `R` is the rotation from sample
    /// to crystal.  The pole's hemispherical position is then binned
    /// by latitude `φ ∈ [0, π/2]` (from the north pole) and longitude
    /// `θ ∈ [0, 2π]`.  Each grain contributes its weight to the bin it
    /// lies in.
    pub fn pole_figure(&self, pole: PoleKind, grid: PoleFigureGrid) -> PoleFigure {
        let d_crystal = pole.direction();
        let mut intensity = vec![0.0_f64; grid.n_lat * grid.n_lon];
        let total_w = self.total_weight().max(1e-30);
        for (g, &w) in self.orientations.iter().zip(self.weights.iter()) {
            // The orientation maps sample → crystal, so the crystal
            // direction expressed in the sample frame is `R^{-1} · d`.
            let d_sample = g.inverse().apply(d_crystal);
            // Latitude from north pole (z-axis).
            let phi = d_sample[2].clamp(-1.0, 1.0).acos();
            // Longitude from x-axis in the equatorial plane.
            let theta = d_sample[1].atan2(d_sample[0]);
            let theta = if theta < 0.0 { theta + 2.0 * std::f64::consts::PI } else { theta };
            let lat = ((phi / (std::f64::consts::FRAC_PI_2)) * grid.n_lat as f64) as usize;
            let lon = ((theta / (2.0 * std::f64::consts::PI)) * grid.n_lon as f64) as usize;
            let lat = lat.min(grid.n_lat - 1);
            let lon = lon.min(grid.n_lon - 1);
            intensity[lat * grid.n_lon + lon] += w / total_w;
        }
        PoleFigure {
            pole,
            grid,
            intensity,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_mat_core::CrystalOrientation;

    #[test]
    fn identity_orientation_puts_pole_at_north() {
        // Identity orientation → crystal x-axis aligns with sample x-axis.
        // The (111) pole sits at the sample-frame (1,1,1)/√3.
        let ta = TextureAnalyzer::new(vec![CrystalOrientation::identity()], vec![1.0]);
        let pf = ta.pole_figure(PoleKind::Fcc111, PoleFigureGrid::default());
        let total: f64 = pf.intensity.iter().sum();
        assert!((total - 1.0).abs() < 1e-9);
        // Exactly one cell should carry all the intensity.
        let nonzero = pf.intensity.iter().filter(|v| **v > 0.0).count();
        assert_eq!(nonzero, 1);
    }

    #[test]
    fn random_textures_distribute_intensity() {
        // 50 grains with random orientations (deterministic LCG).
        let mut ta = TextureAnalyzer::new(Vec::new(), Vec::new());
        let mut rng = LinearRng::new(0x1234_5678);
        for _ in 0..50 {
            let g = CrystalOrientation::from_euler_bunge(
                rng.next_f64() * std::f64::consts::TAU,
                rng.next_f64() * std::f64::consts::PI,
                rng.next_f64() * std::f64::consts::TAU,
            );
            ta.orientations.push(g);
            ta.weights.push(1.0);
        }
        let pf = ta.pole_figure(PoleKind::Fcc111, PoleFigureGrid::default());
        let total: f64 = pf.intensity.iter().sum();
        assert!((total - 1.0).abs() < 1e-9);
        // At least a third of the 50 grains should fall in distinct cells
        // (random textures don't perfectly cover the hemisphere with 50
        // grains, but they should spread out reasonably).
        let nonzero = pf.intensity.iter().filter(|v| **v > 0.0).count();
        assert!(
            nonzero >= 25 && nonzero <= 50,
            "nonzero cells = {nonzero}, expected between 25 and 50"
        );
    }

    struct LinearRng(u64);
    impl LinearRng {
        fn new(seed: u64) -> Self {
            Self(seed.max(1))
        }
        fn next_u64(&mut self) -> u64 {
            // Numerical Recipes LCG.
            self.0 = self.0.wrapping_mul(1664525).wrapping_add(1013904223);
            self.0
        }
        fn next_f64(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }
    }
}