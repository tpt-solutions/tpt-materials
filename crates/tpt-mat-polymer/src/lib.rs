//! Polymer constitutive models.
//!
//! - [`ChainModel`] — freely-jointed chain, worm-like chain,
//!   Arruda–Boyce 8-chain.
//! - [`PolymerModel`] — combined response.
//! - [`stress_strain`] — Cauchy stress under uniaxial extension.

#![warn(missing_docs)]

mod arruda_boyce;
mod chain;
mod fjc;
mod fjc_chain;
mod wlc;

pub use arruda_boyce::{arruda_boyce_nominal_stress, arruda_boyce_true_stress};
pub use chain::{ChainModel, PolymerModel};
pub use fjc::{inverse_langevin_approx, inverse_langevin_exact};
pub use fjc_chain::fjc_force_extension;
pub use wlc::wlc_force_extension;

/// Compute the Cauchy (true) stress for a polymer under
/// uniaxial extension at stretch `λ`.
pub fn stress_strain(model: &PolymerModel, stretch: f64) -> f64 {
    if stretch <= 1.0 {
        return 0.0;
    }
    match &model.chain {
        ChainModel::ArrudaBoyce {
            n_segments,
            shear_modulus,
        } => arruda_boyce_true_stress(*n_segments, *shear_modulus, stretch),
        ChainModel::WormLikeChain {
            persistence_length,
            contour_length,
            k_b,
        } => {
            // Convert the WLC force–extension to a Cauchy stress
            // for a single chain and assume a chain-density of
            // `k_b / k_B T` per unit volume.
            let f = wlc_force_extension(*persistence_length, *contour_length, stretch);
            f * k_b / stretch
        }
        ChainModel::FreelyJointedChain {
            segment_length,
            num_segments,
            k_b,
        } => {
            let f = fjc_force_extension(*segment_length, *num_segments, stretch);
            f * k_b / stretch
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arruda_boyce_stress_grows_with_stretch() {
        use crate::ChainModel;
        let model = PolymerModel {
            chain: ChainModel::ArrudaBoyce {
                n_segments: 8,
                shear_modulus: 0.4e6,
            },
            crosslink_density: 0.0,
        };
        assert!(stress_strain(&model, 2.0) < stress_strain(&model, 4.0));
    }
}
