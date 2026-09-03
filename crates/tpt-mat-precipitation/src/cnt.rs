//! Classical nucleation theory (Turnbull–Fisher steady-state
//! nucleation rate).

/// Boltzmann constant (J/K).
const KB: f64 = 1.380_649e-23;

/// Geometric prefactor `16 π / 3` for the spherical-cap nucleation
/// barrier.
const GEOM_PREFACTOR: f64 = 16.0 * core::f64::consts::PI / 3.0;

/// Spherical-cap homogeneous-nucleation barrier
/// `ΔG* = (16π / 3) γ³ / Δg_v²`.
///
/// All quantities SI:
///
/// - `interfacial_energy` — `γ` (J/m²)
/// - `driving_force_per_volume` — `Δg_v` (J/m³)
pub fn nucleation_barrier(interfacial_energy: f64, driving_force_per_volume: f64) -> f64 {
    let g3 = interfacial_energy.powi(3);
    let dgv2 = driving_force_per_volume.powi(2);
    if dgv2.abs() < 1.0e-30 {
        return f64::INFINITY;
    }
    GEOM_PREFACTOR * g3 / dgv2
}

/// Critical radius `r* = 2 γ / Δg_v` for a spherical nucleus.
pub fn critical_radius(interfacial_energy: f64, driving_force_per_volume: f64) -> f64 {
    let denom = driving_force_per_volume;
    if denom.abs() < 1.0e-30 {
        return f64::INFINITY;
    }
    2.0 * interfacial_energy / denom
}

/// Turnbull–Fisher steady-state nucleation rate
/// `I = I_0 exp(−ΔG*/kT) exp(−Q/kT)`.
///
/// `I_0` is the kinetic prefactor (m⁻³·s⁻¹), `Q` the activation
/// energy for atomic attachment (J), `t_k` the temperature (K),
/// `interfacial_energy` `γ` (J/m²), and `driving_force_per_volume`
/// `Δg_v` (J/m³).
pub fn classical_nucleation_rate(
    i_0: f64,
    activation_energy: f64,
    t_k: f64,
    interfacial_energy: f64,
    driving_force_per_volume: f64,
) -> f64 {
    let barrier = nucleation_barrier(interfacial_energy, driving_force_per_volume);
    let exp_thermo = (-barrier / (KB * t_k)).exp();
    let exp_kinetic = (-activation_energy / (KB * t_k)).exp();
    i_0 * exp_thermo * exp_kinetic
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn critical_radius_decreases_with_driving_force() {
        let r1 = critical_radius(0.5, 1.0e8);
        let r2 = critical_radius(0.5, 1.0e9);
        assert!(r2 < r1);
    }

    #[test]
    fn nucleation_rate_increases_with_driving_force() {
        let i1 = classical_nucleation_rate(1.0e40, 1.0e-19, 600.0, 0.5, 1.0e8);
        let i2 = classical_nucleation_rate(1.0e40, 1.0e-19, 600.0, 0.5, 1.0e9);
        assert!(i2 > i1);
    }

    #[test]
    fn nucleation_rate_grows_with_driving_force_at_constant_temperature() {
        // Use tiny driving forces so the barrier term
        // exp(-ΔG*/kT) doesn't underflow.
        let i_low = classical_nucleation_rate(1.0e40, 1.0e-19, 600.0, 0.5, 1.0e9);
        let i_high = classical_nucleation_rate(1.0e40, 1.0e-19, 600.0, 0.5, 2.0e9);
        assert!(i_high > i_low);
    }

    #[test]
    fn barrier_is_zero_at_zero_interfacial_energy() {
        assert!(approx(nucleation_barrier(0.0, 1.0e8), 0.0, 1.0e-30));
    }
}