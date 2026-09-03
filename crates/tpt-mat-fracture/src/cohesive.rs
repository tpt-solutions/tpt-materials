//! Cohesive-zone traction–separation models.

use serde::{Deserialize, Serialize};

/// Bilinear traction–separation cohesive law.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CohesiveZoneModel {
    /// Peak traction `t_0` (Pa).
    pub peak_traction: f64,
    /// Critical separation `δ_c` (m).
    pub critical_separation: f64,
    /// Fracture energy `G_c` (J/m²).
    pub fracture_energy: f64,
}

impl CohesiveZoneModel {
    /// Construct from `t_0` and `G_c`.  The bilinear law fixes
    /// `δ_c = 2 G_c / t_0` so the area under the curve equals
    /// `G_c`.
    pub fn from_peak_and_fracture_energy(peak_traction: f64, fracture_energy: f64) -> Self {
        let critical_separation = 2.0 * fracture_energy / peak_traction;
        Self {
            peak_traction,
            critical_separation,
            fracture_energy,
        }
    }

    /// Bilinear traction–separation `t(δ)` for `0 ≤ δ ≤ δ_c`.
    /// After `δ_c`, traction is 0 (complete decohesion).
    pub fn traction(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        if delta >= self.critical_separation {
            return 0.0;
        }
        // Linear ramp-up from 0 to t_0 at δ = δ_c / 2, then
        // linear ramp-down from t_0 to 0 at δ = δ_c.
        if delta <= 0.5 * self.critical_separation {
            self.peak_traction * (delta / (0.5 * self.critical_separation))
        } else {
            self.peak_traction
                * ((self.critical_separation - delta)
                    / (0.5 * self.critical_separation))
        }
    }
}

/// A sample of traction vs. separation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TractionSeparation {
    /// Separation `δ` (m).
    pub delta: f64,
    /// Traction `t` (Pa).
    pub traction: f64,
}

/// Benzeggagh–Kenane mixed-mode fracture criterion:
///
/// `G_T = G_I + G_II`, with `G_c = G_Ic + (G_IIc − G_Ic)(G_II/G_T)^η`
/// (Benzeggagh & Kenane 1996).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MixedModeDecomposition {
    /// Mode-I toughness `G_Ic` (J/m²).
    pub g_i_c: f64,
    /// Mode-II toughness `G_IIc` (J/m²).
    pub g_ii_c: f64,
    /// BK exponent `η` (typical 1.5–2.3).
    pub eta: f64,
}

impl MixedModeDecomposition {
    /// Effective mixed-mode toughness `G_c(G_I, G_II)`.
    pub fn effective_toughness(&self, g_i: f64, g_ii: f64) -> f64 {
        let g_i = g_i.max(0.0);
        let g_ii = g_ii.max(0.0);
        let g_total = g_i + g_ii;
        if g_total <= 0.0 {
            return self.g_i_c;
        }
        let ratio = (g_ii / g_total).max(0.0).min(1.0);
        let pow = ratio.powf(self.eta);
        self.g_i_c + (self.g_ii_c - self.g_i_c) * pow
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn traction_zero_below_zero_delta() {
        let c = CohesiveZoneModel::from_peak_and_fracture_energy(1.0e6, 100.0);
        assert_eq!(c.traction(-0.1), 0.0);
    }

    #[test]
    fn traction_zero_above_critical_separation() {
        let c = CohesiveZoneModel::from_peak_and_fracture_energy(1.0e6, 100.0);
        // δ_c = 2 * 100 / 1e6 = 2e-4
        assert_eq!(c.traction(3.0e-4), 0.0);
    }

    #[test]
    fn traction_is_continuous_at_peak() {
        let c = CohesiveZoneModel::from_peak_and_fracture_energy(1.0e6, 100.0);
        let delta_peak = 1.0e-4;
        assert!(approx(c.traction(delta_peak), 1.0e6, 1.0e-6));
    }

    #[test]
    fn cohesive_law_area_equals_g_c() {
        // Integrate the bilinear law numerically with Simpson's
        // rule over many points and confirm we recover G_c.
        let g_c = 200.0;
        let c = CohesiveZoneModel::from_peak_and_fracture_energy(1.0e6, g_c);
        let n = 1001;
        let mut integral = 0.0;
        let mut prev = 0.0;
        for i in 0..n {
            let delta = (i as f64) * c.critical_separation / ((n - 1) as f64);
            let t = c.traction(delta);
            let d_delta = c.critical_separation / ((n - 1) as f64);
            if i == 0 {
                continue;
            }
            integral += 0.5 * (prev + t) * d_delta;
            prev = t;
        }
        assert!(approx(integral, g_c, 1.0e-3));
    }

    #[test]
    fn bk_decomposition_returns_g_i_c_in_pure_mode_i() {
        let bk = MixedModeDecomposition {
            g_i_c: 100.0,
            g_ii_c: 300.0,
            eta: 2.0,
        };
        assert!(approx(bk.effective_toughness(50.0, 0.0), 100.0, 1.0e-9));
    }

    #[test]
    fn bk_decomposition_returns_g_ii_c_in_pure_mode_ii() {
        let bk = MixedModeDecomposition {
            g_i_c: 100.0,
            g_ii_c: 300.0,
            eta: 2.0,
        };
        assert!(approx(bk.effective_toughness(0.0, 50.0), 300.0, 1.0e-9));
    }
}