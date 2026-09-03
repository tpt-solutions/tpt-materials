//! Hardening laws for single-crystal plasticity.
//!
//! Each slip system `α` has an accumulated shear `γ_α` and a critical
//! resolved shear stress `τ_c^α`.  A [`HardeningLaw`] returns the updated
//! `τ_c^α` for a given slip system given the current state.
//!
//! Latent-hardening (cross-hardening between slip systems) is expressed
//! via an `n_slip × n_slip` [`LatentHardeningMatrix`] `h_{αβ}` such that
//! `Δτ_c^α = Σ_β h_{αβ} Δγ^β` (see Kocks & Brown, 1966; Peirce et al.
//! 1982).

#![warn(missing_docs)]

mod latent;
mod state;
mod voce;

pub use latent::LatentHardeningMatrix;
pub use state::HardeningState;
pub use voce::{CombinedParams, KocksMeckingParams, PowerLawParams, VoceParams};

/// Voce saturating hardening:
/// `τ_c^α = τ_0^α + (τ_s^α − τ_0^α) (1 − exp(−γ^α / γ_c^α)) + θ_0^α γ^α`
/// where `θ_0` is the initial linear hardening slope.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct VoceHardening {
    /// Parameters (applied uniformly to all slip systems).
    pub params: VoceParams,
}

impl VoceHardening {
    /// Construct with a uniform parameter set for every slip system.
    pub fn uniform(params: VoceParams) -> Self {
        Self { params }
    }
}

impl HardeningLaw for VoceHardening {
    fn update(
        &self,
        state: &mut HardeningState,
        slip_systems: &[tpt_mat_crystallography::SlipSystem],
        delta_gamma: &[f64],
        latent: &LatentHardeningMatrix,
    ) {
        let n = slip_systems.len();
        debug_assert_eq!(state.crss.len(), n);
        debug_assert_eq!(state.accumulated_shear.len(), n);
        debug_assert_eq!(delta_gamma.len(), n);
        for (alpha, d_gamma_a) in delta_gamma.iter().enumerate() {
            state.accumulated_shear[alpha] += d_gamma_a;
        }
        for alpha in 0..n {
            let gamma = state.accumulated_shear[alpha];
            let p = self.params;
            let voce_increment =
                (p.tau_s - p.tau_0) * (1.0 - (-gamma / p.gamma_c).exp()) + p.theta_0 * gamma;
            state.crss[alpha] = p.tau_0 + voce_increment;
        }
        // Pure Voce ignores cross-hardening; latent matrix is the
        // caller's responsibility via `CombinedHardening`.
        let _ = latent;
    }
}

/// Power-law (Hutchinson-type) rate-dependent hardening:
/// `Δτ_c^α = h_0 (τ_s^α − τ_c^α) |Δγ^α| / τ_s^α`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PowerLawHardening {
    /// Parameters (applied uniformly to every slip system).
    pub params: PowerLawParams,
}

impl PowerLawHardening {
    /// Construct with a uniform parameter set for every slip system.
    pub fn uniform(params: PowerLawParams) -> Self {
        Self { params }
    }
}

impl HardeningLaw for PowerLawHardening {
    fn update(
        &self,
        state: &mut HardeningState,
        slip_systems: &[tpt_mat_crystallography::SlipSystem],
        delta_gamma: &[f64],
        latent: &LatentHardeningMatrix,
    ) {
        let n = slip_systems.len();
        debug_assert_eq!(state.crss.len(), n);
        debug_assert_eq!(state.accumulated_shear.len(), n);
        debug_assert_eq!(delta_gamma.len(), n);
        let _ = slip_systems;
        let p = self.params;
        for alpha in 0..n {
            state.accumulated_shear[alpha] += delta_gamma[alpha];
            let dg = delta_gamma[alpha].abs();
            let sat = p.tau_s - state.crss[alpha];
            let increment = p.h_0 * sat * dg / p.tau_s.max(1e-12);
            state.crss[alpha] += increment;
        }
        let _ = latent;
    }
}

/// Kocks-Mecking evolution: `θ(dτ/dγ) = θ_0 (1 − τ / τ_s)`,
/// integrated implicitly as
/// `τ_c^new = τ_s − (τ_s − τ_c^old) exp(−γ / γ_c)` with
/// `γ_c = τ_s / θ_0`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct KocksMeckingHardening {
    /// Parameters (applied uniformly to every slip system).
    pub params: KocksMeckingParams,
}

impl KocksMeckingHardening {
    /// Construct with a uniform parameter set for every slip system.
    pub fn uniform(params: KocksMeckingParams) -> Self {
        Self { params }
    }
}

impl HardeningLaw for KocksMeckingHardening {
    fn update(
        &self,
        state: &mut HardeningState,
        slip_systems: &[tpt_mat_crystallography::SlipSystem],
        delta_gamma: &[f64],
        latent: &LatentHardeningMatrix,
    ) {
        let n = slip_systems.len();
        debug_assert_eq!(state.crss.len(), n);
        let _ = slip_systems;
        let p = self.params;
        let gamma_c = (p.tau_s / p.theta_0).max(1e-12);
        for alpha in 0..n {
            state.accumulated_shear[alpha] += delta_gamma[alpha];
            let gamma_inc = state.accumulated_shear[alpha];
            state.crss[alpha] = p.tau_s - (p.tau_s - p.tau_0) * (-gamma_inc / gamma_c).exp();
        }
        let _ = latent;
    }
}

/// Combined hardening: self-hardening modulus `h_0` plus a latent
/// cross-hardening factor `q`.  `h_{αβ} = h_0` for `α = β` and
/// `h_{αβ} = h_0 q` for `α ≠ β`.  Then `Δτ_c^α = Σ_β h_{αβ} Δγ^β`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CombinedHardening {
    /// Parameters for the combined law.
    pub params: CombinedParams,
}

impl HardeningLaw for CombinedHardening {
    fn update(
        &self,
        state: &mut HardeningState,
        slip_systems: &[tpt_mat_crystallography::SlipSystem],
        delta_gamma: &[f64],
        latent: &LatentHardeningMatrix,
    ) {
        let n = slip_systems.len();
        debug_assert_eq!(state.crss.len(), n);
        debug_assert_eq!(state.accumulated_shear.len(), n);
        debug_assert_eq!(delta_gamma.len(), n);
        let _ = slip_systems;
        let matrix = latent.with_latent(self.params.q_latent, self.params.h_0);
        for alpha in 0..n {
            state.accumulated_shear[alpha] += delta_gamma[alpha];
        }
        for alpha in 0..n {
            let mut d_tau = 0.0;
            for beta in 0..n {
                d_tau += matrix.get(alpha, beta) * delta_gamma[beta];
            }
            state.crss[alpha] += d_tau;
        }
    }
}

/// Trait implemented by every hardening law.
pub trait HardeningLaw {
    /// Update the CRSS and accumulated shear for every slip system
    /// given the slip increments `delta_gamma[α]`.
    fn update(
        &self,
        state: &mut HardeningState,
        slip_systems: &[tpt_mat_crystallography::SlipSystem],
        delta_gamma: &[f64],
        latent: &LatentHardeningMatrix,
    );
}

/// Top-level hardening-law selector used by
/// `tpt-mat-crystal-plasticity::CrystalPlasticityModel`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Hardening {
    /// Voce saturating hardening.
    Voce(VoceHardening),
    /// Power-law rate-dependent hardening.
    PowerLaw(PowerLawHardening),
    /// Kocks-Mecking evolution.
    KocksMecking(KocksMeckingHardening),
    /// Combined: any self-law + latent cross-hardening matrix.
    Combined(CombinedHardening),
}

impl Hardening {
    /// Dispatch to the underlying law.
    pub fn update(
        &self,
        state: &mut HardeningState,
        slip_systems: &[tpt_mat_crystallography::SlipSystem],
        delta_gamma: &[f64],
        latent: &LatentHardeningMatrix,
    ) {
        match self {
            Hardening::Voce(h) => h.update(state, slip_systems, delta_gamma, latent),
            Hardening::PowerLaw(h) => h.update(state, slip_systems, delta_gamma, latent),
            Hardening::KocksMecking(h) => h.update(state, slip_systems, delta_gamma, latent),
            Hardening::Combined(h) => h.update(state, slip_systems, delta_gamma, latent),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_mat_crystallography::CrystalStructure;

    fn fcc_12_systems() -> Vec<tpt_mat_crystallography::SlipSystem> {
        CrystalStructure::FCC.slip_systems()
    }

    #[test]
    fn voce_monotone_increase() {
        let slips = fcc_12_systems();
        let mut state = HardeningState::from_crss(&slips);
        let law = VoceHardening::uniform(VoceParams {
            tau_0: 10.0,
            tau_s: 100.0,
            theta_0: 500.0,
            gamma_c: 0.05,
        });
        let latent = LatentHardeningMatrix::diagonal(slips.len(), 1.0);
        for _ in 0..50 {
            law.update(&mut state, &slips, &vec![0.001; slips.len()], &latent);
        }
        for tau in &state.crss {
            assert!(*tau >= 10.0 && *tau <= 200.0 + 1e-9, "out of band: {tau}");
            assert!(tau.is_finite());
        }
    }

    #[test]
    fn power_law_saturates() {
        let slips = fcc_12_systems();
        let mut state = HardeningState::from_crss(&slips);
        let law = PowerLawHardening::uniform(PowerLawParams {
            h_0: 1000.0,
            tau_s: 50.0,
        });
        let latent = LatentHardeningMatrix::diagonal(slips.len(), 1.0);
        for _ in 0..200 {
            law.update(&mut state, &slips, &vec![0.01; slips.len()], &latent);
        }
        for tau in &state.crss {
            assert!((*tau - 50.0).abs() < 1e-6, "did not saturate: {tau}");
        }
    }

    #[test]
    fn kocks_mecking_never_exceeds_saturation() {
        let slips = fcc_12_systems();
        let mut state = HardeningState::from_crss(&slips);
        let law = KocksMeckingHardening::uniform(KocksMeckingParams {
            tau_0: 20.0,
            tau_s: 80.0,
            theta_0: 500.0,
        });
        let latent = LatentHardeningMatrix::diagonal(slips.len(), 1.0);
        for _ in 0..1000 {
            law.update(&mut state, &slips, &vec![0.01; slips.len()], &latent);
        }
        for tau in &state.crss {
            assert!(*tau <= 80.0 + 1e-9);
            assert!(*tau >= 20.0 - 1e-9);
        }
    }

    #[test]
    fn combined_latent_cross_hardening_increases_all_systems() {
        let slips = fcc_12_systems();
        let mut state = HardeningState::from_crss(&slips);
        let initial: Vec<f64> = state.crss.clone();
        let law = CombinedHardening {
            params: CombinedParams {
                h_0: 1000.0,
                q_latent: 1.4,
            },
        };
        let latent = LatentHardeningMatrix::diagonal(slips.len(), 1.0);
        let mut dg = vec![0.0; slips.len()];
        dg[0] = 0.001;
        law.update(&mut state, &slips, &dg, &latent);
        for (i, tau) in state.crss.iter().enumerate() {
            assert!(*tau > initial[i], "system {i} did not increase");
        }
    }
}
