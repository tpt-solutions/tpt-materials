//! Hydrogen transport with trapping.

use serde::{Deserialize, Serialize};

/// Oriani local-equilibrium parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OrianiParams {
    /// Lattice diffusivity `D_L` (m²/s).
    pub d_lattice: f64,
    /// Trap density `N_T` (m⁻³).
    pub trap_density: f64,
    /// Binding energy `E_B` (J/mol).
    pub binding_energy_j_per_mol: f64,
    /// Temperature `T` (K).
    pub temperature_k: f64,
}

const R_GAS: f64 = 8.314_462;

/// Effective diffusivity under Oriani local equilibrium:
/// `D_eff = D_L / (1 + ∂C_T/∂C_L)`.
///
/// For a single reversible trap with density `N_T`, occupancy
/// `θ_T` and equilibrium constant `K = exp(−E_B/RT)`:
///
/// `∂C_T/∂C_L = N_T K / (1 + K C_L)²`
pub fn hydrogen_diffusivity_effective(
    c_lattice: f64,
    params: &OrianiParams,
) -> f64 {
    let k_eq = (-params.binding_energy_j_per_mol / (R_GAS * params.temperature_k)).exp();
    let d_trap_dc_l = params.trap_density * k_eq / (1.0 + k_eq * c_lattice).powi(2);
    params.d_lattice / (1.0 + d_trap_dc_l)
}

/// McNabb–Foster kinetic trapping parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct McNabbFosterParams {
    /// Trap density `N_T` (m⁻³).
    pub trap_density: f64,
    /// Forward trapping rate constant `k` (m³/(mol·s)).
    pub k_trap: f64,
    /// Detrapping rate constant `p` (s⁻¹).
    pub p_detrap: f64,
}

/// McNabb–Foster kinetic trapping rate
/// `∂C_T/∂t = k C_L (N_T − C_T) − p C_T`.
///
/// Returns the time derivative of trapped concentration
/// given the current lattice concentration `c_lattice` and
/// trapped concentration `c_trapped`.
pub fn mcnabb_foster_trapping_rate(
    c_lattice: f64,
    c_trapped: f64,
    params: &McNabbFosterParams,
) -> f64 {
    let trap_capacity = params.trap_density - c_trapped;
    let forward = params.k_trap * c_lattice * trap_capacity;
    let backward = params.p_detrap * c_trapped;
    forward - backward
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn trapping_retards_effective_diffusivity() {
        let p = OrianiParams {
            d_lattice: 1.0e-8,
            trap_density: 1.0e25,
            binding_energy_j_per_mol: 30_000.0,
            temperature_k: 300.0,
        };
        let d_eff = hydrogen_diffusivity_effective(1.0e-3, &p);
        assert!(d_eff < p.d_lattice);
    }

    #[test]
    fn zero_traps_yields_d_lattice() {
        let p = OrianiParams {
            d_lattice: 1.0e-8,
            trap_density: 0.0,
            binding_energy_j_per_mol: 30_000.0,
            temperature_k: 300.0,
        };
        let d_eff = hydrogen_diffusivity_effective(1.0e-3, &p);
        assert!(approx(d_eff, p.d_lattice, 1.0e-15));
    }

    #[test]
    fn trapping_rate_zero_at_equilibrium() {
        let p = McNabbFosterParams {
            trap_density: 1.0e25,
            k_trap: 1.0e-3,
            p_detrap: 1.0e-3,
        };
        // If c_L * (N_T - c_T) == c_T (i.e. K = 1), then rate = 0.
        let n_t = p.trap_density;
        let c_t = n_t * p.k_trap / (p.k_trap + p.p_detrap);
        let c_l = c_t * p.p_detrap / (p.k_trap * (n_t - c_t));
        let rate = mcnabb_foster_trapping_rate(c_l, c_t, &p);
        assert!(approx(rate, 0.0, 1.0e-15));
    }

    #[test]
    fn trapping_rate_positive_when_super_saturated() {
        let p = McNabbFosterParams {
            trap_density: 1.0e25,
            k_trap: 1.0e-3,
            p_detrap: 1.0e-3,
        };
        let rate = mcnabb_foster_trapping_rate(1.0e20, 0.0, &p);
        assert!(rate > 0.0);
    }
}