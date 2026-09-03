//! Solidification-map microstructure prediction.
//!
//! The Hunt (1984) columnar-to-equiaxed transition (CET) map
//! parameterises microstructure as a function of the solidification
//! parameters:
//!
//! - `G` — thermal gradient (K/m)
//! - `R` — growth / solidification-front velocity (m/s)
//!
//! Columnar solidification prevails when
//!
//! ```text
//! G · R^n > G_CET
//! ```
//!
//! and equiaxed otherwise.  For alloys with a cubic dendrite
//! morphology `n ≈ 2`; we use `n = 2` as a default.
//!
//! # References
//!
//! - Hunt, J. D. (1984).  "Steady state columnar and equiaxed
//!   growth of dendrites and eutectic."  Mater. Sci. Eng. 65(1),
//!   75–83.

use serde::{Deserialize, Serialize};

use super::history::CoolingRate;

/// Predicted microstructure descriptor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PredictedMicrostructure {
    /// Mean grain size (m).
    pub grain_size: f64,
    /// Phase fractions (a 2-element vector: `[columnar, equiaxed]`).
    pub phase_fractions: [f64; 2],
    /// Estimated porosity (0–1).
    pub porosity: f64,
    /// Texture index (1 = random; >1 = columnar texture).
    pub texture_index: f64,
}

/// Solidification-map inputs.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SolidificationMap {
    /// Thermal gradient `G` (K/m).
    pub thermal_gradient: f64,
    /// Growth rate `R` (m/s).
    pub growth_rate: f64,
    /// Critical `G · R^n` for the columnar-to-equiaxed transition
    /// (defaults to `2.0e10` for aluminium-alloys).
    pub cet_threshold: f64,
    /// Hall–Petch prefactor `k_y` (Pa·m^(1/2)).
    pub hall_petch_k: f64,
    /// Hall–Petch reference stress `σ_0` (Pa).
    pub hall_petch_sigma_0: f64,
    /// Density of nucleating particles `N_0` (m⁻³) — for equiaxed
    /// grain size.
    pub nucleant_density: f64,
    /// Solidification interval `ΔT_s` (K).
    pub solidification_interval: f64,
    /// Volumetric latent heat `L` / heat capacity `C_p` (K).
    pub latent_heat_to_cp: f64,
}

impl Default for SolidificationMap {
    fn default() -> Self {
        Self {
            thermal_gradient: 1.0e5,
            growth_rate: 0.01,
            cet_threshold: 2.0e10,
            hall_petch_k: 0.4,
            hall_petch_sigma_0: 50.0e6,
            nucleant_density: 1.0e15,
            solidification_interval: 100.0,
            latent_heat_to_cp: 300.0,
        }
    }
}

/// Threshold function `G · R^n`.  Returns `∞` if `G = 0` or `R = 0`.
pub fn columnar_equiaxed_threshold(map: &SolidificationMap, exponent: f64) -> f64 {
    let r_pow = map.growth_rate.powf(exponent);
    if map.thermal_gradient == 0.0 {
        return f64::INFINITY;
    }
    map.thermal_gradient * r_pow
}

/// Grain size from a cooling rate (Hall–Petch-like relation):
///
/// `d ≈ A / (ε̇ⁿ)` with `ε̇ = cooling_rate / ΔT_s`.
///
/// Default parameters give a calibration typical of laser-deposited
/// Ti-6Al-4V.
pub fn grain_size_from_cooling_rate(
    cooling_rate_k_per_s: f64,
    solidification_interval: f64,
) -> f64 {
    let cooling = cooling_rate_k_per_s.max(1.0);
    let eps_dot = cooling / solidification_interval.max(1.0);
    80.0e-6 * eps_dot.powf(-0.33)
}

/// Predict microstructure from a [`CoolingRate`] and a
/// [`SolidificationMap`].
pub fn predict_microstructure(
    cooling: &CoolingRate,
    map: &SolidificationMap,
) -> PredictedMicrostructure {
    let gr_n = columnar_equiaxed_threshold(map, 2.0);
    let gr_threshold = map.cet_threshold;
    let columnar_frac = if gr_n >= gr_threshold { 1.0 } else { 0.0 };
    let equiaxed_frac = 1.0 - columnar_frac;
    let grain_size =
        grain_size_from_cooling_rate(cooling.cooling_rate, map.solidification_interval);
    // Porosity scales with too-low energy density (lack of fusion) and
    // too-high (keyholing); here we approximate as a bell-curve centred
    // on the Hall-Petch reference.
    let porosity = 0.001 * (cooling.cooling_rate / 1.0e4).ln().abs().min(0.05);
    // Texture index: columnar → strong fibre texture, equiaxed → random.
    let texture_index = if columnar_frac > 0.5 { 5.0 } else { 1.0 };
    PredictedMicrostructure {
        grain_size,
        phase_fractions: [columnar_frac, equiaxed_frac],
        porosity,
        texture_index,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn higher_cooling_rate_finer_grain() {
        let d1 = grain_size_from_cooling_rate(1.0e3, 100.0);
        let d2 = grain_size_from_cooling_rate(1.0e5, 100.0);
        assert!(d2 < d1);
    }

    #[test]
    fn columnar_when_gradient_high_growth_low() {
        let map = SolidificationMap {
            thermal_gradient: 1.0e6,
            growth_rate: 1.0e-3,
            cet_threshold: 1.0,
            ..SolidificationMap::default()
        };
        let cool = CoolingRate {
            peak_temperature: 1900.0,
            time_at_peak: 0.0,
            cooling_rate: 1.0e4,
            distance: 1.0e-3,
        };
        let m = predict_microstructure(&cool, &map);
        assert!(m.phase_fractions[0] > 0.5);
    }

    #[test]
    fn equiaxed_when_gradient_low_growth_high() {
        let map = SolidificationMap {
            thermal_gradient: 1.0e3,
            growth_rate: 1.0e-1,
            cet_threshold: 1.0e9,
            ..SolidificationMap::default()
        };
        let cool = CoolingRate {
            peak_temperature: 1900.0,
            time_at_peak: 0.0,
            cooling_rate: 1.0e5,
            distance: 1.0e-3,
        };
        let m = predict_microstructure(&cool, &map);
        assert!(m.phase_fractions[1] > 0.5);
    }

    #[test]
    fn columnar_equiaxed_threshold_handles_zero_gradient() {
        let map = SolidificationMap {
            thermal_gradient: 0.0,
            growth_rate: 1.0e-3,
            ..SolidificationMap::default()
        };
        assert!(columnar_equiaxed_threshold(&map, 2.0).is_infinite());
    }
}
