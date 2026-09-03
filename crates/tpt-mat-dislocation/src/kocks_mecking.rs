//! Kocks–Mecking dislocation-density evolution and the Taylor
//! relation between density and flow stress.

use serde::{Deserialize, Serialize};

use super::state::DislocationDensityState;

/// Kocks–Mecking evolution parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct KocksMeckingParams {
    /// Storage coefficient `k₁` (m⁻¹).
    pub k1: f64,
    /// Dynamic-recovery coefficient `k₂` (dimensionless).
    pub k2: f64,
    /// Taylor coefficient `α` (typically 0.3–0.5).
    pub alpha: f64,
    /// Shear modulus `μ` (Pa).
    pub shear_modulus: f64,
    /// Magnitude of the Burgers vector `b` (m).
    pub burgers_vector: f64,
}

/// Forward-Euler Kocks–Mecking step
/// `ρ_{n+1} = ρ_n + (k₁ √ρ_f − k₂ ρ_n) Δγ`,
/// clamped to `ρ ≥ 0`.
///
/// `rho_forest` is the forest density used by the storage term.
pub fn kocks_mecking_step(
    state: &DislocationDensityState,
    deps: f64,
    params: &KocksMeckingParams,
) -> DislocationDensityState {
    let rho = state.rho_ssd.max(0.0);
    let rho_f = state.rho_forest.max(0.0);
    let drho = (params.k1 * rho_f.sqrt() - params.k2 * rho) * deps;
    let rho_new = (rho + drho).max(0.0);
    DislocationDensityState {
        rho_ssd: rho_new,
        rho_gnd: state.rho_gnd,
        rho_forest: rho_new,
    }
}

/// Taylor flow stress `τ = α μ b √ρ` for a single slip system.
pub fn taylor_stress(state: &DislocationDensityState, params: &KocksMeckingParams) -> f64 {
    let rho = state.total().max(0.0);
    params.alpha * params.shear_modulus * params.burgers_vector * rho.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn default_params() -> KocksMeckingParams {
        KocksMeckingParams {
            k1: 7.0e7,
            k2: 5.0,
            alpha: 0.3,
            shear_modulus: 80.0e9,
            burgers_vector: 2.5e-10,
        }
    }

    #[test]
    fn rho_stays_non_negative() {
        let s = DislocationDensityState {
            rho_ssd: 1.0,
            rho_gnd: 0.0,
            rho_forest: 1.0,
        };
        let s2 = kocks_mecking_step(&s, 100.0, &default_params());
        assert!(s2.rho_ssd >= 0.0);
    }

    #[test]
    fn rho_grows_initially_under_small_deps() {
        let s = DislocationDensityState::default();
        let s2 = kocks_mecking_step(&s, 0.01, &default_params());
        assert!(s2.rho_ssd > s.rho_ssd);
    }

    #[test]
    fn rho_saturates_under_extended_deps() {
        let mut s = DislocationDensityState::default();
        let p = default_params();
        for _ in 0..5000 {
            s = kocks_mecking_step(&s, 0.001, &p);
        }
        // After many steps, k1 sqrt(rho_f) ≈ k2 * rho  ⇒ rho ≈ (k1/k2)^2.
        let expected_sat = (p.k1 / p.k2).powi(2);
        assert!(
            approx(s.rho_ssd, expected_sat, expected_sat * 0.1),
            "got {} expected {expected_sat}",
            s.rho_ssd
        );
    }

    #[test]
    fn taylor_stress_increases_with_density() {
        let p = default_params();
        let s_lo = DislocationDensityState {
            rho_ssd: 1.0e12,
            rho_gnd: 0.0,
            rho_forest: 1.0e12,
        };
        let s_hi = DislocationDensityState {
            rho_ssd: 1.0e14,
            rho_gnd: 0.0,
            rho_forest: 1.0e14,
        };
        assert!(taylor_stress(&s_hi, &p) > taylor_stress(&s_lo, &p));
    }
}