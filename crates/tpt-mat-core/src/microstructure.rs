//! Top-level material data container.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::{Composition, Grain, GrainId, Phase, PhaseId};

/// Microstructure-level material description.
///
/// Holds an id + name, the phase inventory, the grain list, the volume
/// element (cube edge length for now), and the bulk temperature.  The
/// EBSD loader and solver crates attach their own derived data through
/// `properties`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MaterialMicrostructure {
    /// Stable id within a project.
    pub id: String,
    /// Human-readable name.
    pub name: String,
    /// All phases, keyed by [`PhaseId`].
    pub phases: HashMap<PhaseId, Phase>,
    /// All grains, keyed by [`GrainId`].
    pub grains: HashMap<GrainId, Grain>,
    /// Volume element edge length (metres).  `0.0` means "unspecified".
    #[serde(default)]
    pub volume_element_edge: f64,
    /// Bulk temperature (Kelvin).  `0.0` means "unspecified".
    #[serde(default)]
    pub temperature: f64,
    /// Free-form additional properties.
    #[serde(default)]
    pub properties: serde_json::Value,
}

impl MaterialMicrostructure {
    /// Construct an empty material with no phases or grains.
    pub fn new(id: impl Into<String>, name: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            name: name.into(),
            phases: HashMap::new(),
            grains: HashMap::new(),
            volume_element_edge: 0.0,
            temperature: 0.0,
            properties: serde_json::Value::Null,
        }
    }

    /// Add a phase.  Replaces any existing phase with the same id.
    pub fn add_phase(&mut self, phase: Phase) -> &mut Self {
        self.phases.insert(phase.id, phase);
        self
    }

    /// Add a grain.  Replaces any existing grain with the same id.
    pub fn add_grain(&mut self, grain: Grain) -> &mut Self {
        self.grains.insert(grain.id, grain);
        self
    }

    /// Look up a phase by id.
    pub fn phase(&self, id: PhaseId) -> Option<&Phase> {
        self.phases.get(&id)
    }

    /// Look up a grain by id.
    pub fn grain(&self, id: GrainId) -> Option<&Grain> {
        self.grains.get(&id)
    }

    /// Total volume fraction across all phases.  Should be 1 for a
    /// closed material.
    pub fn total_volume_fraction(&self) -> f64 {
        self.phases.values().map(|p| p.volume_fraction).sum()
    }

    /// Aggregate composition across all phases (volume-weighted,
    /// atom-fraction output by default).
    pub fn aggregate_composition(
        &self,
        target: crate::CompositionBasis,
    ) -> Result<Composition, crate::CompositionError> {
        if self.phases.is_empty() {
            return Err(crate::CompositionError::Empty);
        }
        let mut all_elements: std::collections::BTreeSet<String> =
            std::collections::BTreeSet::new();
        for p in self.phases.values() {
            for (el, _) in p.composition.iter() {
                all_elements.insert(el.to_string());
            }
        }
        let mut mixed = std::collections::BTreeMap::new();
        for el in &all_elements {
            let mut acc = 0.0;
            for p in self.phases.values() {
                if p.composition.fraction(el).is_some() {
                    let in_moles = p.composition.convert(crate::CompositionBasis::Mole)?;
                    let mf = in_moles.fraction(el).unwrap_or(0.0);
                    acc += p.volume_fraction * mf;
                }
            }
            mixed.insert(el.clone(), acc);
        }
        let as_mole = Composition::new(crate::CompositionBasis::Mole, mixed, 1e-10)?;
        as_mole.convert(target)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CrystalOrientation;

    fn fe_phase() -> Phase {
        let mut m = std::collections::BTreeMap::new();
        m.insert("Fe".into(), 1.0);
        let c = Composition::new(crate::CompositionBasis::Atomic, m, 1e-12).unwrap();
        Phase::new(PhaseId(0), "alpha-Fe", c, 1.0)
    }

    #[test]
    fn empty_total_volume_fraction_is_zero() {
        let m = MaterialMicrostructure::new("mat-0", "test");
        assert_eq!(m.total_volume_fraction(), 0.0);
    }

    #[test]
    fn add_and_lookup() {
        let mut m = MaterialMicrostructure::new("mat-0", "test");
        m.add_phase(fe_phase());
        let g = Grain::new(
            GrainId(0),
            PhaseId(0),
            CrystalOrientation::identity(),
            tpt_math_linalg_fixed::Vec3::ZERO,
            1e-5,
        );
        m.add_grain(g);
        assert!(m.phase(PhaseId(0)).is_some());
        assert!(m.grain(GrainId(0)).is_some());
        assert_eq!(m.total_volume_fraction(), 1.0);
    }
}
