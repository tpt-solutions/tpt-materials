//! Per-grain / per-integration-point state.

use serde::{Deserialize, Serialize};

use tpt_math_linalg_fixed::Mat3;
use tpt_mat_hardening::HardeningState;

/// State of one single crystal during a deformation simulation:
/// accumulated plastic deformation, lattice rotation, and hardening.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SingleCrystalState {
    /// Accumulated plastic deformation tensor `F^p` (initially identity).
    pub plastic_deformation: Mat3,
    /// Lattice rotation `R` from sample frame to current crystal frame
    /// (initially identity).
    pub lattice_rotation: Mat3,
    /// Slip-system hardening state (CRSS + accumulated shear).
    pub hardening: HardeningState,
}

impl SingleCrystalState {
    /// Initialise from the CRSS carried by each slip system.
    pub fn from_crss(crss: &[f64]) -> Self {
        Self {
            plastic_deformation: Mat3::IDENTITY,
            lattice_rotation: Mat3::IDENTITY,
            hardening: HardeningState::new(crss.to_vec(), vec![0.0; crss.len()]),
        }
    }

    /// Number of slip systems.
    pub fn n_slip(&self) -> usize {
        self.hardening.len()
    }
}