//! Redlich–Kister polynomial for the excess Gibbs energy of a binary
//! A–B solution.

use serde::{Deserialize, Serialize};

/// Parameters of a Redlich–Kister polynomial.
///
/// `^E G(x_B) = Σ_ν L_ν (x_A − x_B)^ν` for `ν = 0, 1, …, N − 1`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RedlichKisterParams {
    /// Polynomial coefficients `L_0, L_1, …, L_{N-1}` (energy units).
    pub coefficients: Vec<f64>,
}

impl Default for RedlichKisterParams {
    fn default() -> Self {
        // Symmetric A–B excess; zero first term reproduces an ideal
        // binary.
        Self {
            coefficients: vec![0.0, 0.0],
        }
    }
}

/// Redlich–Kister excess-Gibbs evaluator.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RedlichKister {
    /// Parameters.
    pub params: RedlichKisterParams,
}

impl Default for RedlichKister {
    fn default() -> Self {
        Self {
            params: RedlichKisterParams::default(),
        }
    }
}

impl RedlichKister {
    /// Construct from explicit parameters.
    pub fn new(params: RedlichKisterParams) -> Self {
        Self { params }
    }

    /// `^E G(x_B)` for a binary A–B solution with mole fraction `x_B`.
    /// `x_A = 1 − x_B`.
    pub fn excess(&self, x_b: f64) -> f64 {
        let x_a = 1.0 - x_b;
        let diff = x_a - x_b;
        let mut sum = 0.0;
        let mut power = 1.0_f64;
        for l in &self.params.coefficients {
            sum += l * power;
            power *= diff;
        }
        sum
    }

    /// `d^E G / dx_B` (closed-form using `d(diff^ν)/dx_B = -2 ν diff^{ν-1}`).
    pub fn excess_derivative(&self, x_b: f64) -> f64 {
        let x_a = 1.0 - x_b;
        let diff = x_a - x_b;
        let mut sum = 0.0;
        let mut power = 1.0_f64;
        for (nu_minus_1, l) in self.params.coefficients.iter().enumerate().skip(1) {
            // d/dx_B diff^ν = -2 ν diff^{ν-1}; nu = nu_minus_1 + 1.
            let nu = nu_minus_1 as f64;
            sum += -2.0 * nu * l * power;
            power *= diff;
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_coefficients_is_zero() {
        let rk = RedlichKister::default();
        for &x in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            assert_eq!(rk.excess(x), 0.0);
        }
    }

    #[test]
    fn symmetric_is_unaffected_by_swapping_endpoints() {
        // Use only even-power coefficients so (x_A - x_B)^ν is
        // symmetric about x_B = 0.5 for ν = 0, 2, 4, ...
        let rk = RedlichKister::new(RedlichKisterParams {
            coefficients: vec![1.0, 0.0, 3.0, 0.0, 2.0],
        });
        assert_eq!(rk.excess(0.25), rk.excess(0.75));
        assert_eq!(rk.excess(0.1), rk.excess(0.9));
        assert_eq!(rk.excess(0.4), rk.excess(0.6));
    }

    #[test]
    fn derivative_matches_finite_difference() {
        let rk = RedlichKister::new(RedlichKisterParams {
            coefficients: vec![10.0, 5.0, -2.0, 1.0],
        });
        let h = 1e-6;
        for &x in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let num = (rk.excess(x + h) - rk.excess(x - h)) / (2.0 * h);
            let ana = rk.excess_derivative(x);
            assert!((num - ana).abs() < 1e-4, "x={x} num={num} ana={ana}");
        }
    }
}
