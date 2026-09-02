//! Per-grain hardening state: CRSS and accumulated shear.

use serde::{Deserialize, Serialize};

/// Hardening state for one grain: CRSS and accumulated shear per slip
/// system.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HardeningState {
    /// Current CRSS `τ_c^α` per slip system (MPa).
    pub crss: Vec<f64>,
    /// Accumulated shear `γ^α` per slip system.
    pub accumulated_shear: Vec<f64>,
}

impl HardeningState {
    /// Initialise from the CRSS carried by each [`SlipSystem`](tpt_mat_crystallography::SlipSystem).
    pub fn from_crss(slip_systems: &[tpt_mat_crystallography::SlipSystem]) -> Self {
        let crss: Vec<f64> = slip_systems
            .iter()
            .map(|s| s.critical_resolved_shear_stress)
            .collect();
        Self {
            crss,
            accumulated_shear: vec![0.0; slip_systems.len()],
        }
    }

    /// Construct with explicit CRSS and accumulated-shear vectors.
    pub fn new(crss: Vec<f64>, accumulated_shear: Vec<f64>) -> Self {
        debug_assert_eq!(crss.len(), accumulated_shear.len());
        Self {
            crss,
            accumulated_shear,
        }
    }

    /// Number of slip systems.
    pub fn len(&self) -> usize {
        self.crss.len()
    }

    /// Whether the state has no slip systems.
    pub fn is_empty(&self) -> bool {
        self.crss.is_empty()
    }
}