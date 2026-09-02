//! A single grain in a microstructure.

use serde::{Deserialize, Serialize};

use tpt_math_linalg_fixed::Vec3;

use crate::{Phase, PhaseId};

/// Stable, opaque grain identifier within a [`MaterialMicrostructure`](crate::MaterialMicrostructure).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(transparent)]
pub struct GrainId(pub u32);

/// One grain in a microstructure.
///
/// Stores orientation as a [`CrystalOrientation`](crate::CrystalOrientation),
/// the id of the phase it belongs to, and bookkeeping fields that the
/// Phase 2+ solvers consume.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Grain {
    /// Stable id within a microstructure.
    pub id: GrainId,
    /// Phase this grain is an instance of.
    pub phase_id: PhaseId,
    /// Lattice orientation of the grain (sample → crystal frame).
    pub orientation: crate::CrystalOrientation,
    /// Geometric centroid of the grain in sample-frame coordinates.
    pub centroid: Vec3,
    /// Equivalent sphere radius of the grain (units: metres).
    pub equivalent_radius: f64,
    /// Neighbouring grain ids, if known.  Empty for EBSD datasets that
    /// do not yet have a neighbour table.
    #[serde(default)]
    pub neighbors: Vec<GrainId>,
}

impl Grain {
    /// Build a new grain with no neighbours.
    pub fn new(
        id: GrainId,
        phase_id: PhaseId,
        orientation: crate::CrystalOrientation,
        centroid: Vec3,
        equivalent_radius: f64,
    ) -> Self {
        Self {
            id,
            phase_id,
            orientation,
            centroid,
            equivalent_radius,
            neighbors: Vec::new(),
        }
    }

    /// Convenience: look up the [`Phase`] this grain belongs to in a
    /// material.  Returns `None` if the grain references a phase that
    /// is not in `material`.
    pub fn phase<'a>(&self, material: &'a crate::MaterialMicrostructure) -> Option<&'a Phase> {
        material.phase(self.phase_id)
    }
}
