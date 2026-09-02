//! A single constituent [`Phase`] of a material.

use serde::{Deserialize, Serialize};

use crate::{Composition, CrystalOrientation};

/// Stable, opaque phase identifier within a [`MaterialMicrostructure`](crate::MaterialMicrostructure).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(transparent)]
pub struct PhaseId(pub u32);

/// A single phase of a microstructure.
///
/// Carries the [`Composition`] (always required), a reference to a
/// crystal structure description (added in Phase 2), and free-form
/// physical properties.  Phase 1 only ships the data container; the
/// constitutive-model field arrives in Phase 2.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Phase {
    /// Stable id within a microstructure.
    pub id: PhaseId,
    /// Human-readable name (e.g. `"FCC austenite"`).
    pub name: String,
    /// Composition in some basis.
    pub composition: Composition,
    /// Volume fraction within the microstructure (0, 1].
    pub volume_fraction: f64,
    /// Lattice orientation of the phase reference frame relative to the
    /// sample frame (defaults to identity).
    pub orientation: CrystalOrientation,
    /// Free-form additional properties, serialised as JSON.
    #[serde(default)]
    pub properties: serde_json::Value,
}

impl Phase {
    /// Construct a phase with identity orientation and empty
    /// `properties`.
    pub fn new(
        id: PhaseId,
        name: impl Into<String>,
        composition: Composition,
        volume_fraction: f64,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            composition,
            volume_fraction,
            orientation: CrystalOrientation::identity(),
            properties: serde_json::Value::Null,
        }
    }
}
