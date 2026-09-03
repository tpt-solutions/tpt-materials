//! Flory–Rehner equilibrium swelling (simplified form).
//!
//! The classic Flory–Rehner equation for the polymer volume
//! fraction `φ_p` at swelling equilibrium in a non-ionic gel
//! immersed in a small-molecule solvent is
//!
//! `ln(1 − φ_p) + φ_p + χ φ_p² + ν_e V_1 (φ_p^{1/3} − φ_p/2) = 0`
//!
//! where `ν_e` is the number of elastically active chains per
//! unit dry volume of polymer and `V_1` is the molar volume of
//! the solvent.  This implementation accepts `ν_e V_1` directly
//! as the dimensionless coupling constant `c_ν`.

use serde::{Deserialize, Serialize};

/// Flory–Rehner input parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FloryRehnerParams {
    /// Flory–Huggins polymer–solvent interaction parameter `χ`.
    pub chi: f64,
    /// Dimensionless coupling `ν_e V_1` (chains per dry volume
    /// × molar volume of solvent).
    pub coupling: f64,
}

/// Solve the Flory–Rehner equilibrium for the polymer volume
/// fraction `φ_p ∈ (0, 1)` and return the swelling ratio
/// `Q = 1/φ_p`.
pub fn flory_rehner_swelling_ratio(params: &FloryRehnerParams) -> f64 {
    // Scan a coarse grid to find a bracket where the LHS
    // changes sign.
    let n_scan = 200;
    let mut prev_phi = 1.0e-6;
    let mut prev_lhs = lhs_flory_rehner(prev_phi, params);
    let mut lo = prev_phi;
    let mut hi = 0.999;
    for i in 1..=n_scan {
        let phi = prev_phi
            + (0.999 - 1.0e-6) * (i as f64) / (n_scan as f64);
        let lhs = lhs_flory_rehner(phi, params);
        if prev_lhs * lhs < 0.0 {
            lo = prev_phi;
            hi = phi;
            break;
        }
        prev_phi = phi;
        prev_lhs = lhs;
    }
    for _ in 0..200 {
        let phi = 0.5 * (lo + hi);
        let lhs = lhs_flory_rehner(phi, params);
        if lhs > 0.0 {
            hi = phi;
        } else {
            lo = phi;
        }
        if (hi - lo) < 1.0e-12 * hi.max(1.0e-12) {
            break;
        }
    }
    let phi = 0.5 * (lo + hi);
    1.0 / phi
}

fn lhs_flory_rehner(phi: f64, p: &FloryRehnerParams) -> f64 {
    let one_minus = 1.0 - phi;
    if one_minus <= 0.0 {
        return f64::INFINITY;
    }
    let mix = one_minus.ln() + phi + p.chi * phi * phi;
    let elastic = p.coupling * (phi.powf(1.0 / 3.0) - phi / 2.0);
    mix + elastic
}

/// Convenience wrapper that takes `χ` and `ν_e V_1` directly.
pub fn equilibrium_swelling_ratio(chi: f64, coupling: f64) -> f64 {
    let p = FloryRehnerParams { chi, coupling };
    flory_rehner_swelling_ratio(&p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn swelling_ratio_exceeds_one() {
        let q = equilibrium_swelling_ratio(0.4, 0.1);
        assert!(q > 1.0);
    }

    #[test]
    fn swelling_decreases_with_coupling() {
        // Higher crosslink density -> smaller equilibrium Q.
        let q_lo = equilibrium_swelling_ratio(0.4, 0.05);
        let q_hi = equilibrium_swelling_ratio(0.4, 1.0);
        eprintln!("q_lo={q_lo} q_hi={q_hi}");
        assert!(q_lo > q_hi);
    }

    #[test]
    fn swelling_increases_with_better_solvent() {
        // Lower χ (better solvent) -> more swelling.
        let q_good = equilibrium_swelling_ratio(0.2, 0.1);
        let q_bad = equilibrium_swelling_ratio(0.8, 0.1);
        assert!(q_good > q_bad);
    }
}