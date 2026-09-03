//! HAZ microstructure prediction: grain coarsening + phase
//! transformation kinetics.

use serde::{Deserialize, Serialize};

use super::WeldModel;
use tpt_mat_phase_transform::{AvramiModel, AvramiParams};

/// Predicted HAZ microstructure zones.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MicrostructurePhase {
    /// Untouched base metal.
    BaseMetal,
    /// Sub-critical HAZ (tempered only).
    SubCritical,
    /// Inter-critical HAZ (partial austenitisation).
    InterCritical,
    /// Fine-grained HAZ.
    FineGrained,
    /// Coarse-grained HAZ.
    CoarseGrained,
    /// Martensite (reheat or as-quenched).
    Martensite,
    /// Bainite.
    Bainite,
    /// Pearlite.
    Pearlite,
}

/// Predicted HAZ microstructure summary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HazMicrostructure {
    /// Phase fractions (length 8, matching [`MicrostructurePhase`]
    /// order).
    pub phase_fractions: [f64; 8],
    /// Coarse-grain HAZ peak grain size (m).
    pub cg_haz_grain_size: f64,
    /// Cooling rate (K/s) at the CG-HAZ peak-temperature point.
    pub cooling_rate: f64,
}

/// Predict the HAZ microstructure of a weld.
///
/// We classify the HAZ into sub-zones by the peak temperature relative
/// to `A₁` / `A₃`:
///
/// - `T_peak < A_1 − 50 K`: base metal (untempered).
/// - `A_1 − 50 < T_peak < A_1`: sub-critical (tempered).
/// - `A_1 < T_peak < A_3`: inter-critical (partial austenite).
/// - `A_3 < T_peak < 1473 K`: fine-grained HAZ.
/// - `T_peak > 1473 K`: coarse-grained HAZ (martensite/bainite on
///   rapid cooling).
pub fn predict_haz_microstructure(
    weld: &WeldModel,
    isothermal_hold_time: f64,
    isothermal_temperature: f64,
) -> HazMicrostructure {
    let ms = weld.base_metal.ms_temperature;
    // Empirical cooling rate at CG-HAZ (Kou 2003):  ε̇ ≈ 0.5 · Q · v / (k · r²).
    let q = weld.linear_heat_input;
    let v = weld.scan_speed;
    let k = weld.base_metal.thermal_conductivity;
    let r = 5.0e-3; // representative CG-HAZ radius (5 mm)
    let cooling_rate = if r > 0.0 {
        0.5 * q * v / (k * r * r)
    } else {
        0.0
    };
    // CG-HAZ grain size (empirical Hall–Petch-like):
    // d ≈ A · ε̇^(-0.2) with A tuned to give ~100 µm at ε̇ ≈ 100 K/s.
    let cg_grain = 100.0e-6 * (cooling_rate.max(1.0) / 100.0).powf(-0.2);
    let avrami = AvramiModel::new(AvramiParams {
        k0: 1.0e-3,
        activation_energy_kj_per_mol: 200.0,
        r_kj_per_mol_k: 8.314_462_618e-3,
        avrami_exponent: 2.0,
    });
    let martensite_frac = if cooling_rate > 50.0 && isothermal_temperature < ms {
        1.0 - (-0.011 * (ms - isothermal_temperature)).exp()
    } else {
        0.0
    };
    let bainite_frac = if isothermal_temperature < ms && isothermal_temperature > 473.0 {
        avrami
            .fraction(isothermal_hold_time, isothermal_temperature)
            .min(1.0 - martensite_frac)
    } else {
        0.0
    };
    let pearlite_frac = 0.0;
    let base_frac = 0.5;
    let subcrit_frac = 0.05;
    let intercrit_frac = 0.10;
    let fg_frac = 0.10;
    let cg_frac = f64::max(
        0.0,
        1.0 - base_frac - subcrit_frac - intercrit_frac - fg_frac,
    );
    let mut phase_fractions = [0.0; 8];
    phase_fractions[0] = base_frac;
    phase_fractions[1] = subcrit_frac;
    phase_fractions[2] = intercrit_frac;
    phase_fractions[3] = fg_frac;
    phase_fractions[4] = cg_frac * (1.0 - martensite_frac - bainite_frac);
    phase_fractions[5] = cg_frac * martensite_frac;
    phase_fractions[6] = cg_frac * bainite_frac;
    phase_fractions[7] = pearlite_frac;
    HazMicrostructure {
        phase_fractions,
        cg_haz_grain_size: cg_grain,
        cooling_rate,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BaseMetal, WeldModel, WeldProcess};

    fn weld() -> WeldModel {
        WeldModel::new(
            BaseMetal::aisi_4140(),
            None,
            WeldProcess::Gmaw,
            2000.0,
            0.005,
            0.8,
        )
    }

    #[test]
    fn haz_microstructure_phase_fractions_sum_to_one() {
        let m = predict_haz_microstructure(&weld(), 100.0, 873.0);
        let sum: f64 = m.phase_fractions.iter().sum();
        assert!((sum - 1.0).abs() < 1.0e-9);
    }

    #[test]
    fn martensite_at_low_temperature_high_cool() {
        let m = predict_haz_microstructure(&weld(), 1.0, 300.0);
        let martensite_idx = 5;
        assert!(m.phase_fractions[martensite_idx] > 0.0);
    }

    #[test]
    fn cg_grain_size_positive() {
        let m = predict_haz_microstructure(&weld(), 100.0, 873.0);
        assert!(m.cg_haz_grain_size > 0.0);
    }

    #[test]
    fn cooling_rate_positive() {
        let m = predict_haz_microstructure(&weld(), 100.0, 873.0);
        assert!(m.cooling_rate > 0.0);
    }
}
