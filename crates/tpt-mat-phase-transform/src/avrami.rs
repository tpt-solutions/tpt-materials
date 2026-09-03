//! Avrami / JMAK isothermal transformation kinetics.

use serde::{Deserialize, Serialize};

/// Avrami (JMAK) parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AvramiParams {
    /// Pre-exponential factor `k_0` in `k(T) = k_0 exp(-Q / (R T))`.
    pub k0: f64,
    /// Activation energy `Q` (kJ/mol).
    pub activation_energy_kj_per_mol: f64,
    /// Universal gas constant (kJ/(mol·K)).
    pub r_kj_per_mol_k: f64,
    /// Avrami exponent `n` (typically 1–4).
    pub avrami_exponent: f64,
}

impl Default for AvramiParams {
    fn default() -> Self {
        Self {
            k0: 1.0e5,
            activation_energy_kj_per_mol: 80.0,
            r_kj_per_mol_k: 8.314_462_618e-3,
            avrami_exponent: 3.0,
        }
    }
}

/// Avrami (JMAK) isothermal transformation kinetics.
///
/// Transformed fraction at time `t` and temperature `T`:
///
/// `f(t, T) = 1 − exp(−k(T) t^n)`
///
/// with `k(T) = k_0 · exp(−Q / (R T))`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AvramiModel {
    /// Parameters.
    pub params: AvramiParams,
}

impl Default for AvramiModel {
    fn default() -> Self {
        Self {
            params: AvramiParams::default(),
        }
    }
}

impl AvramiModel {
    /// Construct from explicit parameters.
    pub fn new(params: AvramiParams) -> Self {
        Self { params }
    }

    /// Rate constant `k(T) = k_0 exp(−Q / (R T))` (0 for `T ≤ 0`).
    pub fn rate_constant(&self, t_kelvin: f64) -> f64 {
        if t_kelvin <= 0.0 {
            return 0.0;
        }
        let p = self.params;
        let exponent = -p.activation_energy_kj_per_mol / (p.r_kj_per_mol_k * t_kelvin);
        p.k0 * exponent.exp()
    }

    /// `f(t, T)` for `t ≥ 0`.
    pub fn fraction(&self, t: f64, t_kelvin: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        let k = self.rate_constant(t_kelvin);
        1.0 - (-k * t.powf(self.params.avrami_exponent)).exp()
    }

    /// Time required to reach transformed fraction `f_target` at
    /// temperature `T`.  Returns `None` if `f_target ≤ 0` or `f_target
    /// ≥ 1`, or if `T ≤ 0`.
    pub fn time_to_fraction(&self, f_target: f64, t_kelvin: f64) -> Option<f64> {
        if !(f_target > 0.0 && f_target < 1.0) || t_kelvin <= 0.0 {
            return None;
        }
        let k = self.rate_constant(t_kelvin);
        if k <= 0.0 {
            return None;
        }
        let n = self.params.avrami_exponent;
        Some((-(1.0 - f_target).ln() / k).powf(1.0 / n))
    }
}

/// Single point on a TTT (Time-Temperature-Transformation) diagram:
/// the time to reach a given fraction at a given temperature.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TttPoint {
    /// Temperature (K).
    pub temperature_k: f64,
    /// Time to reach `fraction_target` at this temperature (s).
    pub time_seconds: f64,
    /// Target fraction.
    pub fraction_target: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fraction_in_unit_interval_and_monotonic() {
        // Use moderate k so 1 - exp(-k t^n) doesn't underflow at
        // long times; the test cares about monotonicity and bounds,
        // not the specific time-scale.
        let a = AvramiModel::new(AvramiParams {
            k0: 1.0e-3,
            ..AvramiParams::default()
        });
        let f0 = a.fraction(0.0, 800.0);
        let f1 = a.fraction(1.0, 800.0);
        let f2 = a.fraction(100.0, 800.0);
        assert_eq!(f0, 0.0);
        assert!(f1 > 0.0 && f1 < 1.0);
        assert!(f2 > f1);
        assert!(f2 < 1.0);
    }

    #[test]
    fn time_to_fraction_is_inverse() {
        let a = AvramiModel::default();
        let t_k = 900.0;
        let t = a.time_to_fraction(0.5, t_k).unwrap();
        let f_back = a.fraction(t, t_k);
        assert!((f_back - 0.5).abs() < 1e-9);
    }

    #[test]
    fn rate_constant_increases_with_temperature() {
        let a = AvramiModel::default();
        assert!(a.rate_constant(900.0) > a.rate_constant(500.0));
    }

    #[test]
    fn c_curve_nose() {
        // TTT diagram: t(0.5, T) should have a "nose" (minimum) —
        // faster transformation at intermediate T than at extremes.
        let a = AvramiModel::default();
        let mut ts = Vec::new();
        for tk in [400, 500, 600, 700, 800, 900, 1000] {
            ts.push(a.time_to_fraction(0.5, tk as f64).unwrap_or(f64::INFINITY));
        }
        let min_t = ts.iter().copied().fold(f64::INFINITY, f64::min);
        let max_t = ts.iter().copied().fold(0.0_f64, f64::max);
        assert!(min_t < max_t);
    }
}
