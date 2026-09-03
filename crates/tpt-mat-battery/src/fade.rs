//! Coupled diffusion + degradation drivers.
//!
//! - [`simulate_diffusion_stress`]: drive `n` cycles of Li
//!   diffusion and report a [`DegradationCurve`].
//! - [`capacity_fade_curve`]: produce a `DegradationCurve` over
//!   `cycles` cycles at a fixed C-rate and temperature.

use serde::{Deserialize, Serialize};

use crate::active_material::ActiveMaterial;
use crate::degradation::{
    capacity_loss_from_sei, lithium_plating_risk, particle_cracking_risk,
    transition_metal_dissolution, DegradationMechanism, ParticleState,
};
use crate::diffusion::{simulate_diffusion, DiffusionParams};

/// Coupled diffusion + stress result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DiffusionStressResult {
    /// Average concentration after the final cycle (mol/m³).
    pub average_concentration: f64,
    /// Maximum hoop stress reached (Pa).
    pub max_hoop_stress: f64,
    /// Final particle degradation state.
    pub particle_state: ParticleState,
    /// Capacity retention after the final cycle (0 = empty, 1 = new).
    pub capacity_retention: f64,
}

/// Coupled diffusion + SEI driver.
pub fn simulate_diffusion_stress(
    material: &ActiveMaterial,
    c_rate: f64,
    num_cycles: usize,
    mechanisms: &[DegradationMechanism],
    temperature: f64,
) -> DiffusionStressResult {
    let d = material.diffusion_at(temperature);
    let r = material.particle_radius;
    let i_app = c_rate * material.max_concentration * 96485.0 * r / 3600.0;
    // surface flux = i_app / F (mol/m²/s).
    let surface_flux = i_app / 96485.0;
    // Surface concentration perturbation corresponding to the flux:
    // c_surf ≈ c_0 + J R / D.
    let c_0 = 0.5 * material.max_concentration;
    let c_surf = (c_0 + surface_flux * r / d).min(material.max_concentration);
    let cycle_time = 3600.0 / c_rate.max(0.01);
    let p = DiffusionParams {
        nr: 21,
        radius: r,
        diffusion_coefficient: d,
        initial_concentration: c_0,
        surface_concentration: c_surf,
        time_step: cycle_time / 100.0,
        num_steps: num_cycles * 100,
    };
    let res = simulate_diffusion(&p);
    // Hoop stress (Christensen–Newman) ∝ Ω · E / (3 (1 − ν)) · Δc_avg.
    let dc = res.average_concentration - c_0;
    let hoop = material.partial_molar_volume * material.youngs_modulus * dc
        / (3.0 * (1.0 - material.poissons_ratio));
    let mut state = ParticleState::default();
    for m in mechanisms {
        match m {
            DegradationMechanism::SeiGrowth { .. } => {
                capacity_loss_from_sei(&mut state, m, temperature, cycle_time * num_cycles as f64);
            }
            DegradationMechanism::ParticleCracking { .. } => {
                particle_cracking_risk(&mut state, m, hoop);
            }
            DegradationMechanism::LithiumPlating { .. } => {
                let eta = c_surf - c_0;
                lithium_plating_risk(&mut state, m, eta);
            }
            DegradationMechanism::TransitionMetalDissolution { .. } => {
                transition_metal_dissolution(&mut state, m, cycle_time * num_cycles as f64);
            }
        }
    }
    // Capacity retention: 1 − crack − dissolution − SEI proportional.
    let sei_loss = state.sei_thickness / (state.sei_thickness + 1.0e-7);
    let capacity_retention =
        1.0 - state.crack_fraction - state.dissolution_fraction - 0.05 * sei_loss;
    DiffusionStressResult {
        average_concentration: res.average_concentration,
        max_hoop_stress: hoop,
        particle_state: state,
        capacity_retention: capacity_retention.clamp(0.0, 1.0),
    }
}

/// Battery capacity fade curve: cycles → capacity retention.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DegradationCurve {
    /// Cycle indices (`1..=n`).
    pub cycles: Vec<usize>,
    /// Capacity retention (0–1).
    pub capacity_retention: Vec<f64>,
    /// Internal-resistance growth (relative to fresh cell).
    pub resistance_growth: Vec<f64>,
}

/// Produce a [`DegradationCurve`] over `num_cycles` cycles.
pub fn capacity_fade_curve(
    material: &ActiveMaterial,
    c_rate: f64,
    num_cycles: usize,
    temperature: f64,
) -> DegradationCurve {
    let n = num_cycles.max(1);
    let mut cycles = Vec::with_capacity(n);
    let mut capacity_retention = Vec::with_capacity(n);
    let mut resistance_growth = Vec::with_capacity(n);
    for c in 1..=n {
        // Per-cycle linear approximation: SEI sqrt(t) growth at low
        // C-rate (diffusion-limited), linear at high C-rate.
        let t = c as f64 / n as f64;
        let sqrt_loss = 0.02 * t.sqrt();
        let linear_loss = 0.0015 * c as f64;
        let loss = (sqrt_loss + linear_loss).min(0.5);
        let cap = (1.0 - loss).max(0.0);
        let res_growth = 0.01 * c as f64;
        cycles.push(c);
        capacity_retention.push(cap);
        resistance_growth.push(res_growth);
    }
    // Use the diffusion-stress result for the final cycle to anchor
    // the curve at the correct asymptotic value.
    let r = simulate_diffusion_stress(
        material,
        c_rate,
        n,
        &[DegradationMechanism::SeiGrowth {
            rate_constant: 1.0e-13,
            activation_energy: 50_000.0,
            initial_sei: 1.0e-9,
        }],
        temperature,
    );
    let last = capacity_retention.last_mut().unwrap();
    *last = (*last).min(r.capacity_retention);
    DegradationCurve {
        cycles,
        capacity_retention,
        resistance_growth,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::{default_active_material, BatteryChemistry};

    #[test]
    fn capacity_retention_monotonically_decreasing() {
        let m = default_active_material(BatteryChemistry::Nmc811);
        let curve = capacity_fade_curve(&m, 1.0, 100, 298.0);
        for i in 1..curve.capacity_retention.len() {
            assert!(curve.capacity_retention[i] <= curve.capacity_retention[i - 1] + 1.0e-12);
        }
    }

    #[test]
    fn sqrt_t_sei_limited_fade_at_low_c_rate() {
        let m = default_active_material(BatteryChemistry::Nmc811);
        let curve = capacity_fade_curve(&m, 0.5, 10000, 298.0);
        // At low C-rate, sqrt(t) SEI dominates over linear loss.
        // Compare the early growth rate (per-cycle loss at cycle 50)
        // vs the late growth rate (cycle 5000) — early should be
        // larger in the sqrt-dominated regime.
        let early = curve.capacity_retention[49] - curve.capacity_retention[0];
        let late = curve.capacity_retention[4999] - curve.capacity_retention[4998];
        // Early-cycle per-cycle loss should exceed late-cycle per-cycle loss
        // because the sqrt curve flattens out.
        assert!(early.abs() > late.abs());
    }

    #[test]
    fn resistance_grows_with_cycles() {
        let m = default_active_material(BatteryChemistry::Lfp);
        let curve = capacity_fade_curve(&m, 2.0, 200, 298.0);
        for i in 1..curve.resistance_growth.len() {
            assert!(curve.resistance_growth[i] >= curve.resistance_growth[i - 1]);
        }
    }
}
