//! Koistinen–Marburger diffusionless (martensitic) kinetics.

use serde::{Deserialize, Serialize};

/// Koistinen–Marburger parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct KMParams {
    /// Martensite-start temperature `M_s` (K).
    pub m_start_k: f64,
    /// Koistinen–Marburger exponent `α` (1/K).
    /// Typical value for steels: `α ≈ 0.011` (K^-1).
    pub alpha_per_k: f64,
    /// As-quenched retained austenite fraction at `M_s` (1 - f(M_s))
    /// — typically a small positive number that accounts for the
    /// finite kinetics at the start temperature.
    pub retained_at_m_start: f64,
}

impl Default for KMParams {
    fn default() -> Self {
        Self {
            m_start_k: 600.0,
            alpha_per_k: 0.011,
            retained_at_m_start: 0.01,
        }
    }
}

/// Koistinen–Marburger transformation kinetics.
///
/// `f(T) = 1 − exp(−α (M_s − T))` for `T ≤ M_s`, and `f = 0` for
/// `T > M_s`.  Optionally scales the retained fraction at `M_s` to
/// honour a non-zero `retained_at_m_start`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct KoistinenMarburger {
    /// Parameters.
    pub params: KMParams,
}

impl Default for KoistinenMarburger {
    fn default() -> Self {
        Self {
            params: KMParams::default(),
        }
    }
}

impl KoistinenMarburger {
    /// Construct from explicit parameters.
    pub fn new(params: KMParams) -> Self {
        Self { params }
    }

    /// Transformed fraction at temperature `T`.
    pub fn fraction(&self, t_kelvin: f64) -> f64 {
        let p = self.params;
        if t_kelvin >= p.m_start_k {
            return 0.0;
        }
        let f = 1.0 - (-p.alpha_per_k * (p.m_start_k - t_kelvin)).exp();
        // Cap at the maximum recoverable fraction (1 - retained at M_s).
        let max_f = 1.0 - p.retained_at_m_start;
        f.min(max_f).clamp(0.0, 1.0)
    }

    /// Derivative `df/dT` at `T` (negative below `M_s`, 0 above).
    pub fn derivative(&self, t_kelvin: f64) -> f64 {
        let p = self.params;
        if t_kelvin >= p.m_start_k {
            return 0.0;
        }
        // d(1 - exp(-α(M_s - T)))/dT = -α · exp(-α(M_s - T))
        -p.alpha_per_k * (-p.alpha_per_k * (p.m_start_k - t_kelvin)).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_transformation_above_ms() {
        let km = KoistinenMarburger::default();
        assert_eq!(km.fraction(km.params.m_start_k + 1.0), 0.0);
        assert_eq!(km.fraction(km.params.m_start_k + 100.0), 0.0);
    }

    #[test]
    fn fraction_increases_on_cooling() {
        let km = KoistinenMarburger::default();
        let t1 = km.params.m_start_k - 10.0;
        let t2 = km.params.m_start_k - 100.0;
        let t3 = km.params.m_start_k - 200.0;
        let f1 = km.fraction(t1);
        let f2 = km.fraction(t2);
        let f3 = km.fraction(t3);
        assert!(f2 > f1);
        assert!(f3 > f2);
        assert!(f3 < 1.0);
    }

    #[test]
    fn retained_at_ms_honoured() {
        // retained_at_m_start caps the maximum recoverable fraction:
        // the fraction should never exceed 1 - retained.
        let km = KoistinenMarburger {
            params: KMParams {
                m_start_k: 600.0,
                alpha_per_k: 0.5,
                retained_at_m_start: 0.1,
            },
        };
        let f_low_t = km.fraction(km.params.m_start_k - 100.0);
        assert!(f_low_t <= 1.0 - km.params.retained_at_m_start + 1e-12);
    }

    #[test]
    fn derivative_sign() {
        let km = KoistinenMarburger::default();
        assert!(km.derivative(km.params.m_start_k - 50.0) < 0.0);
        assert_eq!(km.derivative(km.params.m_start_k + 1.0), 0.0);
    }
}
