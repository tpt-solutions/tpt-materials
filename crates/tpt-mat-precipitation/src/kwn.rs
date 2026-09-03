//! Kampmann–Wagner numerical (KWN) size-class solver.
//!
//! Discretise the precipitate size distribution into `N`
//! classes with radii `R[i]`.  Forward-Euler update:
//!
//! ```text
//! n_i(t + Δt) = n_i(t) + (J_{i-1/2} − J_{i+1/2}) Δt + I_nuc δ(i − i_seed) Δt
//! ```
//!
//! where the flux `J_{i+1/2} = n_i (dR/dt)_i` is the upwind
//! flux evaluated at the class interfaces.  This implementation
//! uses a simple explicit forward-Euler scheme — sufficient for
//! the calibration tests in the suite; production runs should
//! switch to a fully-implicit or operator-split scheme.

use serde::{Deserialize, Serialize};

/// KWN solver parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KwnParams {
    /// Discrete class radii `R[i]` (m), monotonically increasing.
    pub class_radii: Vec<f64>,
    /// Bulk supersaturation `S = (c − c_eq) / c_eq` driving the
    /// growth rate `dR/dt = K_g (S − R*/R)` (simplified Zener).
    pub supersaturation: f64,
    /// Kinetic constant `K_g` (m/s).
    pub growth_constant: f64,
    /// Capillary length `R*` (m), typically a few nm.
    pub capillary_radius: f64,
}

impl KwnParams {
    /// Construct from radii and Zener coefficients.
    pub fn new(
        class_radii: Vec<f64>,
        supersaturation: f64,
        growth_constant: f64,
        capillary_radius: f64,
    ) -> Self {
        Self {
            class_radii,
            supersaturation,
            growth_constant,
            capillary_radius,
        }
    }

    /// Per-class growth rate `dR/dt = K_g (S − R*/R)`.
    pub fn growth_rate(&self, r: f64) -> f64 {
        if r <= 0.0 {
            return 0.0;
        }
        self.growth_constant * (self.supersaturation - self.capillary_radius / r)
    }
}

/// KWN state: number density per class (m⁻³).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KwnState {
    /// Number density per size class (m⁻³).
    pub number_densities: Vec<f64>,
}

impl KwnState {
    /// Empty state matching a given `class_radii` length.
    pub fn empty_for(params: &KwnParams) -> Self {
        Self {
            number_densities: vec![0.0; params.class_radii.len()],
        }
    }

    /// Total number of precipitates (m⁻³).
    pub fn total_number(&self) -> f64 {
        self.number_densities.iter().sum()
    }

    /// Volume fraction (assuming spherical precipitates and
    /// unit volume fraction per particle `= (4/3) π R³`).
    pub fn volume_fraction(&self, params: &KwnParams) -> f64 {
        let mut vf = 0.0;
        for (n, r) in self
            .number_densities
            .iter()
            .zip(params.class_radii.iter())
        {
            vf += n * (4.0 / 3.0) * core::f64::consts::PI * r.powi(3);
        }
        vf
    }
}

/// One explicit KWN time step.  `nucleation_flux` (m⁻³·s⁻¹)
/// seeds the smallest class.
pub fn kwn_step(
    state: &KwnState,
    dt: f64,
    params: &KwnParams,
    nucleation_flux: f64,
) -> KwnState {
    let n_classes = params.class_radii.len();
    if n_classes < 2 {
        return state.clone();
    }
    let mut next = state.number_densities.clone();
    // Upwind fluxes: J_{i+1/2} = n_i * dR/dt at class i.
    let mut flux = vec![0.0; n_classes + 1];
    for i in 0..n_classes {
        let r = params.class_radii[i];
        flux[i + 1] = state.number_densities[i] * params.growth_rate(r);
    }
    for i in 0..n_classes {
        next[i] += (flux[i] - flux[i + 1]) * dt;
        next[i] = next[i].max(0.0);
    }
    if nucleation_flux > 0.0 {
        next[0] += nucleation_flux * dt;
    }
    KwnState {
        number_densities: next,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_params() -> KwnParams {
        KwnParams::new(
            (1..=10).map(|i| (i as f64) * 1.0e-9).collect(),
            0.05,
            1.0e-9,
            2.0e-9,
        )
    }

    #[test]
    fn empty_state_has_zero_total() {
        let s = KwnState::empty_for(&test_params());
        assert_eq!(s.total_number(), 0.0);
    }

    #[test]
    fn zero_dt_returns_same_state() {
        let mut s = KwnState::empty_for(&test_params());
        s.number_densities[5] = 1.0e20;
        let s2 = kwn_step(&s, 0.0, &test_params(), 0.0);
        assert_eq!(s, s2);
    }

    #[test]
    fn nucleation_flux_seeds_smallest_class() {
        let s = KwnState::empty_for(&test_params());
        let s2 = kwn_step(&s, 1.0, &test_params(), 1.0e20);
        assert!(s2.number_densities[0] > 0.0);
    }

    #[test]
    fn volume_fraction_grows_with_seed() {
        let s = KwnState::empty_for(&test_params());
        let s2 = kwn_step(&s, 1.0, &test_params(), 1.0e20);
        let p = test_params();
        assert!(s2.volume_fraction(&p) > 0.0);
    }
}