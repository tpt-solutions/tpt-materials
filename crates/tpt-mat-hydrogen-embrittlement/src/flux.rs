//! Stress-driven hydrogen flux.
//!
//! The chemical potential of hydrogen in a stressed solid has
//! a hydrostatic-stress contribution `-V_H σ_h`, leading to an
//! extra uphill-diffusion flux `J = -D C V_H / (RT) ∇σ_h` in
//! addition to the Fickian `J = -D ∇C` term.

const R_GAS: f64 = 8.314_462;

/// Compute the stress-driven flux density `J_stress = -D C V_H / (RT) * ∂σ_h/∂x`.
///
/// # Arguments
///
/// - `diffusivity` — `D` (m²/s)
/// - `concentration` — `C` (mol/m³)
/// - `partial_molar_volume` — `V_H` (m³/mol)
/// - `temperature_k` — `T` (K)
/// - `gradient_hydrostatic_stress` — `∂σ_h/∂x` (Pa/m)
pub fn stress_driven_flux(
    diffusivity: f64,
    concentration: f64,
    partial_molar_volume: f64,
    temperature_k: f64,
    gradient_hydrostatic_stress: f64,
) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    -diffusivity * concentration * partial_molar_volume / (R_GAS * temperature_k)
        * gradient_hydrostatic_stress
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn zero_gradient_yields_zero_flux() {
        let f = stress_driven_flux(1.0e-8, 1.0, 2.0e-6, 300.0, 0.0);
        assert!(approx(f, 0.0, 1.0e-30));
    }

    #[test]
    fn positive_gradient_yields_negative_flux() {
        // H accumulates at regions of high hydrostatic tension.
        let f = stress_driven_flux(1.0e-8, 1.0, 2.0e-6, 300.0, 1.0e3);
        assert!(f < 0.0);
    }

    #[test]
    fn zero_temperature_yields_zero_flux() {
        let f = stress_driven_flux(1.0e-8, 1.0, 2.0e-6, 0.0, 1.0e3);
        assert!(approx(f, 0.0, 1.0e-30));
    }
}