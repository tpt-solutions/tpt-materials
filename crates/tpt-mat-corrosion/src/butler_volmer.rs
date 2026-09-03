//! Butler–Volmer charge-transfer kinetics and Tafel asymptotes.
//!
//! The net current density across a single electrochemical
//! half-reaction is
//!
//! ```text
//! i(E) = i_0 · ( exp( (E − E_eq) / b_a ) − exp( −(E − E_eq) / b_c ) )
//! ```
//!
//! where `b_a = RT / (α_a F n)` and `b_c = RT / (α_c F n)` are the
//! anodic and cathodic Tafel slopes (V/dec, where "dec" means
//! decade, i.e. the natural log).  Here `R = 8.314 J/mol·K`,
//! `F = 96485 C/mol`, `α_a, α_c` are the charge-transfer
//! coefficients (`α_a + α_c = 1`), and `n` is the number of
//! electrons.
//!
//! # References
//!
//! - Butler, J. A. V. (1924).  "Studies in heterogeneous
//!   equilibria. Part II.—The kinetic interpretation of the
//!   Nernst theory of electromotive force."  Trans. Faraday Soc.
//!   19, 729–733.
//! - Erdey-Grúz, T., & Volmer, M. (1930).  "Zur Theorie der
//!   Wasserstoffüberspannung."  Z. Phys. Chem. 150A, 203–213.

use serde::{Deserialize, Serialize};

use crate::{FARADAY, GAS_CONSTANT};

/// Compute the Tafel slope `b = RT / (α n F)` in V/decade.
///
/// Inputs:
/// - `temperature` — T (K)
/// - `alpha` — charge-transfer coefficient (`0 < α < 1`)
/// - `n` — number of electrons transferred.
pub fn tafel_slope(temperature: f64, alpha: f64, n: f64) -> f64 {
    (GAS_CONSTANT * temperature) / (alpha * n * FARADAY)
}

/// Tafel anodic current density: `i_a = i_0 exp(η / b_a)`.
///
/// `η = E − E_eq` is the overpotential.
pub fn tafel_anodic(i_0: f64, b_a: f64, overpotential: f64) -> f64 {
    i_0 * (overpotential / b_a).exp()
}

/// Tafel cathodic current density magnitude:
/// `|i_c| = i_0 exp(−η / b_c)`.
pub fn tafel_cathodic(i_0: f64, b_c: f64, overpotential: f64) -> f64 {
    i_0 * (-overpotential / b_c).exp()
}

/// A single electrode half-reaction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ElectrodeKinetics {
    /// Equilibrium potential `E_eq` (V vs SHE).
    pub equilibrium_potential: f64,
    /// Exchange current density `i_0` (A/m²).
    pub exchange_current_density: f64,
    /// Anodic Tafel slope `b_a` (V/dec, natural log).
    pub tafel_slope_anodic: f64,
    /// Cathodic Tafel slope `b_c` (V/dec, natural log).
    pub tafel_slope_cathodic: f64,
}

impl ElectrodeKinetics {
    /// Construct from temperature, charge-transfer coefficients and
    /// `n` (number of electrons).  Convenience helper that builds the
    /// Tafel slopes via `RT / (α F n)`.
    pub fn from_alphas(
        equilibrium_potential: f64,
        exchange_current_density: f64,
        alpha_a: f64,
        alpha_c: f64,
        n: f64,
        temperature: f64,
    ) -> Self {
        Self {
            equilibrium_potential,
            exchange_current_density,
            tafel_slope_anodic: tafel_slope(temperature, alpha_a, n),
            tafel_slope_cathodic: tafel_slope(temperature, alpha_c, n),
        }
    }
}

/// Butler–Volmer current density (A/m²) for an `ElectrodeKinetics`.
///
/// `i(E) = i_0 [ exp(η / b_a) − exp(−η / b_c) ]`,
/// `η = E − E_eq`.
pub fn butler_volmer_current_density(kinetics: &ElectrodeKinetics, potential: f64) -> f64 {
    let eta = potential - kinetics.equilibrium_potential;
    butler_volmer_current(
        kinetics.exchange_current_density,
        kinetics.tafel_slope_anodic,
        kinetics.tafel_slope_cathodic,
        eta,
    )
}

/// Butler–Volmer current with explicit Tafel slopes.
pub fn butler_volmer_current(i_0: f64, b_a: f64, b_c: f64, overpotential: f64) -> f64 {
    i_0 * ((overpotential / b_a).exp() - (-overpotential / b_c).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn tafel_slope_rt_alpha_nF() {
        // α=0.5, n=1, T=298 K ⇒ b = 8.314·298 / (0.5·1·96485) = 0.0514 V/dec.
        let b = tafel_slope(298.0, 0.5, 1.0);
        assert!(approx(b, 0.0514, 1.0e-3));
    }

    #[test]
    fn tafel_anodic_grows_exponentially() {
        let i = tafel_anodic(1.0e-3, 0.05, 0.10);
        assert!(i > 1.0e-3);
        let i2 = tafel_anodic(1.0e-3, 0.05, 0.20);
        assert!(i2 > i);
    }

    #[test]
    fn butler_volmer_zero_at_equilibrium() {
        let k = ElectrodeKinetics {
            equilibrium_potential: 0.0,
            exchange_current_density: 1.0e-4,
            tafel_slope_anodic: 0.05,
            tafel_slope_cathodic: 0.05,
        };
        assert!(approx(butler_volmer_current_density(&k, 0.0), 0.0, 1.0e-15));
    }

    #[test]
    fn butler_volmer_anodic_positive() {
        let k = ElectrodeKinetics {
            equilibrium_potential: 0.0,
            exchange_current_density: 1.0e-4,
            tafel_slope_anodic: 0.05,
            tafel_slope_cathodic: 0.05,
        };
        assert!(butler_volmer_current_density(&k, 0.10) > 0.0);
    }

    #[test]
    fn butler_volmer_cathodic_negative() {
        let k = ElectrodeKinetics {
            equilibrium_potential: 0.0,
            exchange_current_density: 1.0e-4,
            tafel_slope_anodic: 0.05,
            tafel_slope_cathodic: 0.05,
        };
        assert!(butler_volmer_current_density(&k, -0.10) < 0.0);
    }
}
