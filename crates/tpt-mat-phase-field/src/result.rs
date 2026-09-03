//! Phase-field result snapshot.

use serde::{Deserialize, Serialize};

/// Result of a phase-field step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhaseFieldResult {
    /// Order-parameter field(s).  Allen-Cahn / Kobayashi: length 1.
    /// Cahn-Hilliard: length 1 (concentration).  Multi-phase:
    /// `num_grains` flat arrays of length `n_cells`.
    pub order_parameter: Vec<Vec<f64>>,
    /// Concentration field (Cahn-Hilliard only).
    #[serde(default)]
    pub concentration: Vec<f64>,
    /// Dimensionless temperature field (Kobayashi only).
    #[serde(default)]
    pub temperature: Vec<f64>,
    /// Total free energy at this snapshot.
    pub free_energy: f64,
    /// Interface area (count of cells where `|∇η| > threshold`).
    pub interface_area: f64,
    /// Time reached at this snapshot.
    pub time: f64,
}
