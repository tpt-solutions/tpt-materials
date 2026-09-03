//! Data-source provenance per spec §9.

use serde::{Deserialize, Serialize};

/// Provenance for a material property value.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "value")]
pub enum DataSource {
    /// Peer-reviewed or book reference.
    Textbook(String),
    /// ASTM / ISO standard.
    Standard(String),
    /// NIST database.
    Nist(String),
    /// Internal lab measurement (with measurement ID).
    LabMeasurement {
        /// Measurement identifier.
        measurement_id: String,
    },
    /// Manufacturer datasheet.
    Datasheet(String),
}