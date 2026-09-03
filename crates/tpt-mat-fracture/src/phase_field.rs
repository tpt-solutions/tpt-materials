//! Phase-field fracture: AT1 / AT2 degradation functions and
//! crack-surface dissipation density.

use serde::{Deserialize, Serialize};

/// Phase-field fracture model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PhaseFieldFracture {
    /// Bourdin–Francfort–Marigo / AT2 model: `w(d) = d²`,
    /// `g(d) = (1 − d)² + κ` (small residual stiffness `κ`).
    At2,
    /// Ambati–Amor–Bourdin / AT1 model: `w(d) = d`,
    /// `g(d) = (1 − d)²`.
    At1,
}

/// AT1 degradation function `g(d) = (1 − d)²`.
///
/// `g(0) = 1` (undamaged), `g(1) = 0` (fully broken).
pub fn at1_degradation(d: f64) -> f64 {
    let x = 1.0 - d;
    x * x
}

/// AT2 degradation function `g(d) = (1 − d)² + κ` with the
/// residual stiffness `κ` (a small positive number).
pub fn at2_degradation(d: f64) -> f64 {
    at1_degradation(d) + KAPPA_RESIDUAL
}

/// Residual stiffness used by [`at2_degradation`].
const KAPPA_RESIDUAL: f64 = 1.0e-10;

/// Crack-surface dissipation density
/// `(G_c / c_w) [w(d) + l_0² |∇d|²]` per unit volume.
///
/// - `g_c`  — fracture energy (J/m²)
/// - `c_w`  — normalisation factor (`8/3` for AT2, `2` for AT1)
/// - `l_0`  — regularisation length (m)
/// - `d`    — phase-field damage variable (`0 ≤ d ≤ 1`)
/// - `grad_d_sq` — squared gradient `|∇d|²` (1/m²)
pub fn dissipation_density(
    model: PhaseFieldFracture,
    g_c: f64,
    l_0: f64,
    d: f64,
    grad_d_sq: f64,
) -> f64 {
    let c_w = match model {
        PhaseFieldFracture::At1 => 2.0,
        PhaseFieldFracture::At2 => 8.0 / 3.0,
    };
    let w = match model {
        PhaseFieldFracture::At1 => d,
        PhaseFieldFracture::At2 => d * d,
    };
    let l_0_sq = l_0 * l_0;
    (g_c / c_w) * (w + l_0_sq * grad_d_sq)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn at1_degradation_is_one_at_d_zero() {
        assert!(approx(at1_degradation(0.0), 1.0, 1.0e-12));
    }

    #[test]
    fn at1_degradation_is_zero_at_d_one() {
        assert!(approx(at1_degradation(1.0), 0.0, 1.0e-12));
    }

    #[test]
    fn at1_degradation_monotonically_decreases() {
        let mut last = at1_degradation(0.0);
        for i in 1..=10 {
            let d = i as f64 / 10.0;
            let v = at1_degradation(d);
            assert!(v < last);
            last = v;
        }
    }

    #[test]
    fn dissipation_is_zero_for_undamaged_uniform_field() {
        // d = 0, |∇d|² = 0 ⇒ dissipation = 0.
        let d = dissipation_density(PhaseFieldFracture::At2, 100.0, 0.01, 0.0, 0.0);
        assert!(approx(d, 0.0, 1.0e-12));
    }

    #[test]
    fn dissipation_is_maximum_at_full_damage() {
        let d = dissipation_density(PhaseFieldFracture::At1, 100.0, 0.01, 1.0, 0.0);
        // (G_c / 2) * 1 = 50
        assert!(approx(d, 50.0, 1.0e-9));
    }

    #[test]
    fn dissipation_includes_gradient_term() {
        let d0 = dissipation_density(PhaseFieldFracture::At2, 100.0, 0.01, 0.5, 0.0);
        let d1 = dissipation_density(PhaseFieldFracture::At2, 100.0, 0.01, 0.5, 100.0);
        let l_0_sq = 0.01 * 0.01;
        let expected_gradient_term = (100.0 / (8.0 / 3.0)) * l_0_sq * 100.0;
        assert!(approx(d1 - d0, expected_gradient_term, 1.0e-6));
    }
}
