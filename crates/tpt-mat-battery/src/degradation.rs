//! Degradation mechanisms for Li-ion battery active materials.
//!
//! Each [`DegradationMechanism`] returns a [`ParticleState`]
//! capturing the *current* loss state of a single particle.  The
//! [`DegradationCurve`] driver in [`crate::fade`] composes the
//! mechanisms over cycles.

use serde::{Deserialize, Serialize};

/// A degradation mechanism.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum DegradationMechanism {
    /// SEI growth on the anode / positive-electrode particle surface
    /// (Peled 1979).  Solved as
    ///
    /// `dδ/dt = k · exp(−Q / RT)`
    ///
    /// with `δ` (m) the SEI thickness; consumed Li scales linearly
    /// with `δ`.
    SeiGrowth {
        /// SEI growth rate constant `k` (m/s).
        rate_constant: f64,
        /// Activation energy `Q` (J/mol).
        activation_energy: f64,
        /// Initial SEI thickness (m).
        initial_sei: f64,
    },
    /// Particle cracking: instant loss of connectivity when the
    /// diffusion-induced stress exceeds the material strength.
    ParticleCracking {
        /// Critical stress `σ_c` (Pa).
        critical_stress: f64,
    },
    /// Lithium plating: parasitic Li deposition when the surface
    /// potential falls below 0 V vs Li/Li⁺.
    LithiumPlating {
        /// Plating potential (V vs Li/Li⁺) below which plating occurs.
        plating_potential: f64,
    },
    /// Transition-metal dissolution: first-order decay of the
    /// active material.
    TransitionMetalDissolution {
        /// Dissolution rate `k` (1/s).
        dissolution_rate: f64,
    },
}

/// Per-particle degradation state.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ParticleState {
    /// SEI thickness `δ` (m).
    pub sei_thickness: f64,
    /// 0 = intact, 1 = fully cracked.
    pub crack_fraction: f64,
    /// Fraction of the active material dissolved (0 = none).
    pub dissolution_fraction: f64,
    /// `true` if lithium plating is currently active.
    pub is_plating: bool,
}

impl Default for ParticleState {
    fn default() -> Self {
        Self {
            sei_thickness: 0.0,
            crack_fraction: 0.0,
            dissolution_fraction: 0.0,
            is_plating: false,
        }
    }
}

const R_GAS: f64 = 8.314_462_618;

/// Update the SEI thickness.
///
/// `δ(t + dt) = δ(t) + k · exp(−Q / RT) · dt`.
pub fn capacity_loss_from_sei(
    sei: &mut ParticleState,
    mechanism: &DegradationMechanism,
    temperature: f64,
    time_increment: f64,
) {
    if let DegradationMechanism::SeiGrowth {
        rate_constant,
        activation_energy,
        initial_sei,
    } = mechanism
    {
        if sei.sei_thickness == 0.0 {
            sei.sei_thickness = *initial_sei;
        }
        let growth =
            rate_constant * (-activation_energy / (R_GAS * temperature)).exp() * time_increment;
        sei.sei_thickness += growth;
    }
}

/// Decide whether the diffusion-induced hoop stress exceeds the
/// critical cracking stress.
pub fn particle_cracking_risk(
    state: &mut ParticleState,
    mechanism: &DegradationMechanism,
    hoop_stress: f64,
) {
    if let DegradationMechanism::ParticleCracking { critical_stress } = mechanism {
        if hoop_stress >= *critical_stress && state.crack_fraction < 1.0 {
            state.crack_fraction = 1.0;
        }
    }
}

/// Detect lithium plating from the surface overpotential.
pub fn lithium_plating_risk(
    state: &mut ParticleState,
    mechanism: &DegradationMechanism,
    surface_overpotential: f64,
) {
    if let DegradationMechanism::LithiumPlating { plating_potential } = mechanism {
        state.is_plating = surface_overpotential < *plating_potential;
    }
}

/// First-order transition-metal dissolution update.
pub fn transition_metal_dissolution(
    state: &mut ParticleState,
    mechanism: &DegradationMechanism,
    time_increment: f64,
) {
    if let DegradationMechanism::TransitionMetalDissolution { dissolution_rate } = mechanism {
        let decay = dissolution_rate * time_increment;
        state.dissolution_fraction = (state.dissolution_fraction + decay).min(1.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sei_grows_arrhenius() {
        let mut s = ParticleState::default();
        let m = DegradationMechanism::SeiGrowth {
            rate_constant: 1.0e-13,
            activation_energy: 50_000.0,
            initial_sei: 1.0e-9,
        };
        capacity_loss_from_sei(&mut s, &m, 298.0, 3600.0);
        assert!(s.sei_thickness > 1.0e-9);
    }

    #[test]
    fn sei_grows_faster_at_higher_temperature() {
        let m = DegradationMechanism::SeiGrowth {
            rate_constant: 1.0e-13,
            activation_energy: 50_000.0,
            initial_sei: 1.0e-9,
        };
        let mut s1 = ParticleState::default();
        capacity_loss_from_sei(&mut s1, &m, 298.0, 3600.0);
        let mut s2 = ParticleState::default();
        capacity_loss_from_sei(&mut s2, &m, 318.0, 3600.0);
        assert!(s2.sei_thickness > s1.sei_thickness);
    }

    #[test]
    fn particle_cracks_when_stress_exceeds_critical() {
        let mut s = ParticleState::default();
        let m = DegradationMechanism::ParticleCracking {
            critical_stress: 500.0e6,
        };
        particle_cracking_risk(&mut s, &m, 600.0e6);
        assert!((s.crack_fraction - 1.0).abs() < 1.0e-9);
        particle_cracking_risk(&mut s, &m, 100.0e6);
        assert!((s.crack_fraction - 1.0).abs() < 1.0e-9);
    }

    #[test]
    fn lithium_plating_at_low_overpotential() {
        let mut s = ParticleState::default();
        let m = DegradationMechanism::LithiumPlating {
            plating_potential: 0.0,
        };
        lithium_plating_risk(&mut s, &m, -0.05);
        assert!(s.is_plating);
        lithium_plating_risk(&mut s, &m, 0.05);
        assert!(!s.is_plating);
    }

    #[test]
    fn tm_dissolution_first_order() {
        let mut s = ParticleState::default();
        let m = DegradationMechanism::TransitionMetalDissolution {
            dissolution_rate: 1.0e-9,
        };
        transition_metal_dissolution(&mut s, &m, 86400.0 * 365.0);
        assert!(s.dissolution_fraction > 0.0);
        assert!(s.dissolution_fraction < 1.0);
    }
}
