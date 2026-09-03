//! Crystallographic-slip (Tanaka–Mura 1981) fatigue-indicator parameter.
//!
//! The Tanaka–Mura dislocation model says a fatigue crack initiates on
//! a slip band once the stored dislocation energy reaches a critical
//! value, equivalent to
//!
//! ```text
//! Σ_α γ_acc^α  ≥  γ_c
//! ```
//!
//! where `γ_acc^α` is the accumulated shear on slip system `α` and
//! `γ_c` is a material constant.  In the Coffin–Manson form the
//! cycles-to-initiation from a single slip system obeys
//!
//! ```text
//! γ_acc^α = Δγ^α · N   ⇒   N_i = γ_c / Δγ^α
//! ```
//!
//! # Reference
//!
//! - Tanaka, K., & Mura, T. (1981).  "A dislocation model for
//!   fatigue crack initiation."  J. Appl. Mech. 48(1), 97–103.

/// `true` iff the accumulated shear on a slip system has reached the
/// critical value for crack initiation.
pub fn critical_shear_initiation(accumulated_shear: f64, critical: f64) -> bool {
    accumulated_shear >= critical
}

/// Cycles-to-initiation from a Coffin–Manson fatigue law on a single
/// slip system:
///
/// ```text
/// γ_acc^α = γ'_f · (2 N_i)^c   ⇒   N_i = (1/2) (γ_acc / γ'_f)^(1/c)
/// ```
///
/// Inputs:
///
/// - `delta_gamma` — plastic shear strain amplitude on the dominant
///   slip system
/// - `gamma_f_prime` — fatigue-ductility coefficient on the slip plane
/// - `c` — fatigue-ductility exponent (negative, typically `-0.6`).
///
/// Returns the cycles-to-initiation; `∞` if `c == 0` or
/// `delta_gamma == 0`.
pub fn cycles_to_initiation_coffin_manson(delta_gamma: f64, gamma_f_prime: f64, c: f64) -> f64 {
    if c == 0.0 || delta_gamma <= 0.0 || gamma_f_prime <= 0.0 {
        return f64::INFINITY;
    }
    // (γ'_f · (2 N)^c) = Δγ ⇒ 2 N = (Δγ / γ'_f)^(1/c)
    let ratio = delta_gamma / gamma_f_prime;
    if ratio <= 0.0 {
        return f64::INFINITY;
    }
    let two_n = ratio.powf(1.0 / c);
    0.5 * two_n
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn critical_accumulated_shear_initiation() {
        let cum = vec![0.05, 0.10, 0.20, 0.18, 0.04];
        let crit = 0.15;
        let mut predicted = Vec::new();
        for &c in &cum {
            predicted.push(critical_shear_initiation(c, crit));
        }
        assert!(!predicted[0]);
        assert!(predicted[2]);
        assert!(predicted[3]);
        assert!(!predicted[4]);
    }

    #[test]
    fn cycles_to_initiation_coffin_manson_returns_finite() {
        let nf = cycles_to_initiation_coffin_manson(0.01, 0.5, -0.6);
        let nf2 = cycles_to_initiation_coffin_manson(0.02, 0.5, -0.6);
        assert!(nf > 0.0);
        assert!(nf > nf2);
    }

    #[test]
    fn cycles_to_initiation_zero_returns_infinite() {
        let nf = cycles_to_initiation_coffin_manson(0.0, 0.5, -0.6);
        assert!(nf.is_infinite());
    }
}
