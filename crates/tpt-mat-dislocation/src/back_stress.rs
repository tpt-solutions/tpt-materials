//! Armstrong–Frederick kinematic back-stress update.

use serde::{Deserialize, Serialize};

/// Armstrong–Frederick kinematic-hardening parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BackStressParams {
    /// Kinematic-hardening modulus `c` (Pa).
    pub c_modulus: f64,
    /// Dynamic-recovery coefficient `γ` (dimensionless).
    pub gamma: f64,
}

/// One forward-Euler Armstrong–Frederick step:
///
/// `X_{n+1} = X_n + c Δε^p − γ X_n |Δε^p|`.
///
/// `delta_eps_p` is the accumulated plastic-strain increment
/// over the step (signed scalar).
pub fn armstrong_frederick_step(
    back_stress: f64,
    delta_eps_p: f64,
    params: &BackStressParams,
) -> f64 {
    back_stress + params.c_modulus * delta_eps_p - params.gamma * back_stress * delta_eps_p.abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn zero_recovery_grows_back_stress() {
        let p = BackStressParams {
            c_modulus: 100.0e6,
            gamma: 0.0,
        };
        let x = armstrong_frederick_step(0.0, 0.01, &p);
        assert!(approx(x, 1.0e6, 1.0e-6));
    }

    #[test]
    fn recovery_saturates_back_stress() {
        // After many steps at large Δε^p, X → c / γ.
        let p = BackStressParams {
            c_modulus: 100.0e6,
            gamma: 100.0,
        };
        let mut x = 0.0;
        for _ in 0..1000 {
            x = armstrong_frederick_step(x, 0.01, &p);
        }
        let expected = p.c_modulus / p.gamma;
        assert!(approx(x, expected, expected * 0.05));
    }

    #[test]
    fn load_reversal_reduces_back_stress() {
        let p = BackStressParams {
            c_modulus: 100.0e6,
            gamma: 100.0,
        };
        let x_forward = armstrong_frederick_step(0.0, 0.01, &p);
        let x_after_reverse = armstrong_frederick_step(x_forward, -0.02, &p);
        assert!(x_after_reverse < x_forward);
    }
}
