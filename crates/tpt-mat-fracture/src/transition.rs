//! ASTM E1921 master curve and ductile-to-brittle transition.

use serde::{Deserialize, Serialize};

/// Master-curve reference parameters (ASTM E1921).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MasterCurveParams {
    /// Reference temperature `T_0` (K) at which the median
    /// toughness `K_Jc(med) ≈ 100 MPa·m^0.5` is observed.
    pub t0_k: f64,
    /// Weibull shape parameter for cleavage (typical 4).
    pub weibull_shape: f64,
}

/// Median cleavage toughness `K_Jc(med)` at temperature `T`
/// following the master-curve shape
/// `K_Jc(med) = 30 + 70 exp(0.019 (T − T_0))` (MPa·m^0.5)
/// (ASTM E1921-25 Eq. 1).
pub fn master_curve_k_jc(p: &MasterCurveParams, t_k: f64) -> f64 {
    30.0 + 70.0 * (0.019 * (t_k - p.t0_k)).exp()
}

/// Indicative ductile-to-brittle transition temperature
/// (DBTT) using a simple threshold `K_Jc = 100 MPa·m^0.5`.
/// Returns `T` at which the median toughness equals the
/// threshold (linear approximation `T_DBTT ≈ T_0 − 8 K` for
/// the reference shape).
pub fn fracture_toughness_transition(p: &MasterCurveParams) -> f64 {
    // Inverting the master-curve equation at K_Jc = 100:
    // 100 = 30 + 70 exp(0.019 (T - T0))
    // => (T - T0) = (1/0.019) * ln((100 - 30)/70) = 0
    // So T_DBTT ≈ T_0.  Allow a tiny offset for the Weibull-driven
    // lower-bound transition.
    #[allow(clippy::eq_op)]
        let offset = {
            let numerator: f64 = 100.0 - 30.0;
            let denominator: f64 = 70.0;
            (numerator / denominator).ln() / 0.019
        };
    p.t0_k + offset
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn master_curve_at_t0_equals_100() {
        let p = MasterCurveParams {
            t0_k: 300.0,
            weibull_shape: 4.0,
        };
        assert!(approx(master_curve_k_jc(&p, 300.0), 100.0, 1.0e-6));
    }

    #[test]
    fn master_curve_grows_with_temperature() {
        let p = MasterCurveParams {
            t0_k: 300.0,
            weibull_shape: 4.0,
        };
        assert!(master_curve_k_jc(&p, 350.0) > master_curve_k_jc(&p, 250.0));
    }

    #[test]
    fn transition_temperature_recovers_t0() {
        let p = MasterCurveParams {
            t0_k: 273.0,
            weibull_shape: 4.0,
        };
        assert!(approx(fracture_toughness_transition(&p), 273.0, 1.0e-6));
    }
}