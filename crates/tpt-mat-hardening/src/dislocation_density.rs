//! Dislocation-density-based hardening law.
//!
//! Implements the Kocks–Mecking evolution
//! `dρ/dγ = k₁ √ρ_f − k₂ ρ` per slip system, with the Taylor
//! relation `τ = α μ b √ρ` providing the CRSS.  State is held in a
//! `tpt-mat-dislocation::DislocationDensityState` per slip system.

use serde::{Deserialize, Serialize};

use tpt_mat_crystallography::SlipSystem;
use tpt_mat_dislocation::{
    kocks_mecking_step, taylor_stress, DislocationDensityState, KocksMeckingParams,
};

use crate::{HardeningLaw, HardeningState, LatentHardeningMatrix};

/// Parameters for the dislocation-density-based hardening law.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DislocationDensityParams {
    /// Kocks–Mecking storage coefficient `k₁` (m⁻¹).
    pub k1: f64,
    /// Kocks–Mecking dynamic-recovery coefficient `k₂` (dimensionless).
    pub k2: f64,
    /// Taylor coefficient `α` (typically 0.3–0.5).
    pub alpha: f64,
    /// Reference shear modulus `μ` (Pa).
    pub shear_modulus: f64,
    /// Magnitude of the Burgers vector `b` (m).
    pub burgers_vector: f64,
    /// Initial CRSS `τ_0` (Pa) at zero dislocation density.
    pub tau_0: f64,
}

impl Default for DislocationDensityParams {
    fn default() -> Self {
        Self {
            k1: 7.0e7,
            k2: 5.0,
            alpha: 0.3,
            shear_modulus: 80.0e9,
            burgers_vector: 2.5e-10,
            tau_0: 30.0e6,
        }
    }
}

/// Dislocation-density-based hardening law.  Drives the CRSS via the
/// Taylor relation from a per-slip-system density state that is
/// updated by the Kocks–Mecking evolution equation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DislocationDensityHardening {
    /// Law parameters.
    pub params: DislocationDensityParams,
}

impl Default for DislocationDensityHardening {
    fn default() -> Self {
        Self {
            params: DislocationDensityParams::default(),
        }
    }
}

impl DislocationDensityHardening {
    /// Construct from explicit parameters.
    pub fn new(params: DislocationDensityParams) -> Self {
        Self { params }
    }

    /// Per-slip-system density state.  Initialised to zero at every
    /// slip system.
    pub fn initial_densities(&self, slips: &[SlipSystem]) -> Vec<DislocationDensityState> {
        vec![DislocationDensityState::default(); slips.len()]
    }

    fn km_params(&self) -> KocksMeckingParams {
        KocksMeckingParams {
            k1: self.params.k1,
            k2: self.params.k2,
            alpha: self.params.alpha,
            shear_modulus: self.params.shear_modulus,
            burgers_vector: self.params.burgers_vector,
        }
    }
}

impl HardeningLaw for DislocationDensityHardening {
    fn update(
        &self,
        state: &mut HardeningState,
        slip_systems: &[SlipSystem],
        delta_gamma: &[f64],
        _latent: &LatentHardeningMatrix,
    ) {
        if slip_systems.is_empty() {
            return;
        }
        let mut densities: Vec<DislocationDensityState> =
            if state.extra.len() == slip_systems.len() * 3 {
                let mut out = Vec::with_capacity(slip_systems.len());
                for chunk in state.extra.chunks_exact(3) {
                    out.push(DislocationDensityState {
                        rho_ssd: chunk[0],
                        rho_gnd: chunk[1],
                        rho_forest: chunk[2],
                    });
                }
                out
            } else {
                self.initial_densities(slip_systems)
            };

        let km = self.km_params();
        for (alpha_idx, slip) in slip_systems.iter().enumerate() {
            let dep = delta_gamma.get(alpha_idx).copied().unwrap_or(0.0).abs();
            if dep > 0.0 {
                densities[alpha_idx] = kocks_mecking_step(&densities[alpha_idx], dep, &km);
            }
            let _ = slip;
            let tau = self.params.tau_0 + taylor_stress(&densities[alpha_idx], &km);
            if let Some(slot) = state.crss.get_mut(alpha_idx) {
                *slot = tau;
            }
        }
        let mut flat = Vec::with_capacity(densities.len() * 3);
        for d in &densities {
            flat.push(d.rho_ssd);
            flat.push(d.rho_gnd);
            flat.push(d.rho_forest);
        }
        state.extra = flat;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::HardeningState;
    use tpt_mat_crystallography::CrystalStructure;

    #[test]
    fn density_law_updates_crs_via_taylor_relation() {
        let slips = CrystalStructure::FCC.slip_systems();
        let mut state = HardeningState::from_crss(&slips);
        let law = DislocationDensityHardening::default();
        let latent = LatentHardeningMatrix::diagonal(slips.len(), 1.0);
        for _ in 0..100 {
            law.update(&mut state, &slips, &vec![0.001; slips.len()], &latent);
        }
        for &tau in &state.crss {
            assert!(tau > law.params.tau_0, "CRSS should exceed tau_0: {tau}");
            assert!(tau.is_finite());
        }
    }

    #[test]
    fn density_law_saturates_when_k1_equals_k2_times_sqrt_rho() {
        let slips = CrystalStructure::FCC.slip_systems();
        let mut state = HardeningState::from_crss(&slips);
        let params = DislocationDensityParams {
            k1: 1.0e8,
            k2: 10.0,
            alpha: 0.3,
            shear_modulus: 80.0e9,
            burgers_vector: 2.5e-10,
            tau_0: 0.0,
        };
        let law = DislocationDensityHardening::new(params);
        let latent = LatentHardeningMatrix::diagonal(slips.len(), 1.0);
        for _ in 0..5000 {
            law.update(&mut state, &slips, &vec![0.001; slips.len()], &latent);
        }
        let expected_rho = (params.k1 / params.k2).powi(2);
        let expected_tau =
            params.alpha * params.shear_modulus * params.burgers_vector * expected_rho.sqrt();
        for &tau in &state.crss {
            let rel = (tau - expected_tau).abs() / expected_tau;
            assert!(rel < 0.1, "saturation CRSS {tau} vs theory {expected_tau}");
        }
    }
}
