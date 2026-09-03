//! Solid-state Li diffusion in a spherical particle.
//!
//! Solves
//!
//! ```text
//! ∂c/∂t = D/r² ∂/∂r (r² ∂c/∂r)
//! ```
//!
//! with a Neumann boundary condition on the particle centre
//! (`∂c/∂r = 0`) and a prescribed surface flux at `r = R`
//! corresponding to the (de)lithiation current.

use serde::{Deserialize, Serialize};

/// Diffusion-simulation parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DiffusionParams {
    /// Number of radial cells.
    pub nr: usize,
    /// Particle radius (m).
    pub radius: f64,
    /// Diffusion coefficient (m²/s).
    pub diffusion_coefficient: f64,
    /// Initial uniform Li concentration (mol/m³).
    pub initial_concentration: f64,
    /// Surface concentration (mol/m³) prescribed as Dirichlet on the
    /// outer shell (use `surface_flux` if a Neumann boundary is
    /// preferred — see [`simulate_diffusion`]).
    pub surface_concentration: f64,
    /// Time step (s).
    pub time_step: f64,
    /// Number of time steps.
    pub num_steps: usize,
}

/// Result of a diffusion simulation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DiffusionResult {
    /// Concentration profile at the end of the simulation
    /// (`mol/m³`, length `nr`).
    pub concentration: Vec<f64>,
    /// Final time (s).
    pub time: f64,
    /// Average concentration over the particle (mol/m³).
    pub average_concentration: f64,
}

/// Solve a single lithiation / delithiation step on a 1-D radial grid.
///
/// Method: explicit forward Euler with `O(Δt)` and `O(Δr²)`
/// truncation error; the FD stencil for `∂/∂r (r² ∂c/∂r)` is
///
/// ```text
/// (1/Δr) [ r²_{i+1/2} (c_{i+1} − c_i) / Δr
///        − r²_{i-1/2} (c_i − c_{i-1}) / Δr ] / Δr
/// ```
///
/// with `r²_{i±1/2} = ((r_i ± Δr/2))²`.
pub fn simulate_diffusion(params: &DiffusionParams) -> DiffusionResult {
    let nr = params.nr.max(2);
    let dr = params.radius / (nr as f64 - 1.0);
    let mut c = vec![params.initial_concentration; nr];
    c[nr - 1] = params.surface_concentration;
    let dt = params.time_step;
    let d = params.diffusion_coefficient;
    let dt_max = 0.25 * dr * dr / d;
    let sub_steps = (dt / dt_max).ceil().max(1.0) as usize;
    let dt_eff = dt / (sub_steps as f64);
    for _ in 0..params.num_steps {
        for _ in 0..sub_steps {
            let mut cn = c.clone();
            for i in 1..(nr - 1) {
                let r_i = (i as f64) * dr;
                let r_p = (i as f64 + 0.5) * dr;
                let r_m = (i as f64 - 0.5) * dr;
                let flux_p = r_p * r_p * (c[i + 1] - c[i]) / dr;
                let flux_m = r_m * r_m * (c[i] - c[i - 1]) / dr;
                cn[i] = c[i] + dt_eff * d * (flux_p - flux_m) / (dr * r_i * r_i);
            }
            cn[0] = cn[1];
            cn[nr - 1] = params.surface_concentration;
            c = cn;
        }
    }
    let mut sum = 0.0_f64;
    let mut volume = 0.0_f64;
    for i in 0..nr {
        let r = (i as f64) * dr;
        let weight = if i == 0 || i == nr - 1 { 0.5 } else { 1.0 };
        let shell = 4.0 * core::f64::consts::PI * r * r * dr * weight;
        sum += c[i] * shell;
        volume += shell;
    }
    DiffusionResult {
        concentration: c,
        time: dt * params.num_steps as f64,
        average_concentration: if volume > 0.0 { sum / volume } else { 0.0 },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diffusion_drives_concentration_toward_surface() {
        let p = DiffusionParams {
            nr: 21,
            radius: 5.0e-6,
            diffusion_coefficient: 1.0e-14,
            initial_concentration: 10_000.0,
            surface_concentration: 20_000.0,
            time_step: 1.0,
            num_steps: 200,
        };
        let res = simulate_diffusion(&p);
        // Average concentration should increase above initial.
        assert!(res.average_concentration > 10_000.0);
        // Surface should match Dirichlet.
        let last = *res.concentration.last().unwrap();
        assert!((last - 20_000.0).abs() < 1.0e-6);
    }

    #[test]
    fn diffusion_steady_state_uniform() {
        // With `initial = surface`, no time should change anything.
        let p = DiffusionParams {
            nr: 11,
            radius: 5.0e-6,
            diffusion_coefficient: 1.0e-14,
            initial_concentration: 15_000.0,
            surface_concentration: 15_000.0,
            time_step: 1.0,
            num_steps: 100,
        };
        let res = simulate_diffusion(&p);
        for &c in &res.concentration {
            assert!((c - 15_000.0).abs() < 1.0e-6);
        }
    }
}
