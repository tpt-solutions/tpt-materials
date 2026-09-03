//! Arrhenius temperature dependence of diffusion coefficients.

use serde::{Deserialize, Serialize};

/// Arrhenius-law diffusivity parameters.
///
/// `D(T) = D_0 · exp(-Q / (R T))`
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ArrheniusParams {
    /// Pre-exponential factor `D_0` (units match the diffusion length²
    /// / time convention used by the caller; here we just expose it as
    /// a dimensionless multiplier — physical units are not enforced).
    pub d0: f64,
    /// Activation energy `Q` in kJ/mol.
    pub activation_energy_kj_per_mol: f64,
    /// Universal gas constant `R` in kJ/(mol·K).  Defaults to
    /// `8.314462618e-3` kJ/(mol·K) so callers can pass `T` in K.
    pub r_kj_per_mol_k: f64,
}

impl Default for ArrheniusParams {
    fn default() -> Self {
        Self {
            d0: 1.0,
            activation_energy_kj_per_mol: 100.0,
            r_kj_per_mol_k: 8.314_462_618e-3,
        }
    }
}

/// Arrhenius-law diffusivity evaluator.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ArrheniusDiffusivity {
    /// Parameters.
    pub params: ArrheniusParams,
}

impl Default for ArrheniusDiffusivity {
    fn default() -> Self {
        Self {
            params: ArrheniusParams::default(),
        }
    }
}

impl ArrheniusDiffusivity {
    /// Construct from explicit parameters.
    pub fn new(params: ArrheniusParams) -> Self {
        Self { params }
    }

    /// `D(T)` (returns 0 for `T ≤ 0`).
    pub fn at_temperature(&self, t_kelvin: f64) -> f64 {
        if t_kelvin <= 0.0 {
            return 0.0;
        }
        let p = self.params;
        let exponent = -p.activation_energy_kj_per_mol / (p.r_kj_per_mol_k * t_kelvin);
        p.d0 * exponent.exp()
    }

    /// Temperature derivative `dD/dT = D · Q / (R T²)` (positive).
    pub fn temperature_derivative(&self, t_kelvin: f64) -> f64 {
        let d = self.at_temperature(t_kelvin);
        if t_kelvin <= 0.0 {
            return 0.0;
        }
        d * self.params.activation_energy_kj_per_mol
            / (self.params.r_kj_per_mol_k * t_kelvin * t_kelvin)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arrhenius_increases_with_temperature() {
        let a = ArrheniusDiffusivity::default();
        let d_low = a.at_temperature(500.0);
        let d_high = a.at_temperature(1000.0);
        assert!(d_high > d_low);
        assert!(d_low > 0.0);
    }

    #[test]
    fn derivative_positive() {
        let a = ArrheniusDiffusivity::default();
        assert!(a.temperature_derivative(800.0) > 0.0);
    }

    #[test]
    fn zero_kelvin_returns_zero() {
        let a = ArrheniusDiffusivity::default();
        assert_eq!(a.at_temperature(0.0), 0.0);
    }
}
