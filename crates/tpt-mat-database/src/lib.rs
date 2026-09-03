//! Materials property database (Phase 8, spec §5 Domain 10).
//!
//! - [`MaterialRecord`]: one row of mechanical / thermal / electrical
//!   properties plus provenance.
//! - [`MaterialsDatabase`]: ordered set of records with JSON I/O.
//! - [`PropertyQuery`]: range query (min/max) over a single property.
//! - [`DataSource`]: traceable provenance (spec §9).
//!
//! The bundled data set is intentionally small and uses values that
//! are either textbook reference data or that the project can
//! publish under MIT/Apache-2.0.  Larger datasets should be loaded
//! from external JSON files via [`MaterialsDatabase::from_json`].

#![warn(missing_docs)]

mod builtin;
mod provenance;
mod query;
mod record;

pub use provenance::DataSource;
pub use query::{PropertyQuery, Property};
pub use record::{
    Composition, Electrical, Mechanical, MaterialRecord, Thermal,
};

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Ordered set of material records.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct MaterialsDatabase {
    /// Records, in insertion order.
    pub materials: Vec<MaterialRecord>,
}

/// Database errors.
#[derive(Debug, Error)]
pub enum DatabaseError {
    /// JSON could not be parsed.
    #[error("invalid JSON: {0}")]
    Json(#[from] serde_json::Error),
}

impl MaterialsDatabase {
    /// Construct an empty database.
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct the bundled MIT-clean database.
    pub fn load_builtin() -> Self {
        Self {
            materials: builtin::builtin_records(),
        }
    }

    /// Add a record.
    pub fn push(&mut self, record: MaterialRecord) {
        self.materials.push(record);
    }

    /// Number of records.
    pub fn len(&self) -> usize {
        self.materials.len()
    }

    /// `true` iff no records are present.
    pub fn is_empty(&self) -> bool {
        self.materials.is_empty()
    }

    /// All materials whose `name` exactly matches.
    pub fn search_by_name(&self, name: &str) -> Vec<&MaterialRecord> {
        self.materials
            .iter()
            .filter(|m| m.name == name)
            .collect()
    }

    /// All materials matching a property-range query.
    pub fn search_by_property(&self, q: PropertyQuery) -> Vec<&MaterialRecord> {
        self.materials
            .iter()
            .filter(|m| q.matches(m))
            .collect()
    }

    /// Parse from JSON (the format used by `to_json`).
    pub fn from_json(s: &str) -> Result<Self, DatabaseError> {
        Ok(serde_json::from_str(s)?)
    }

    /// Serialise to JSON.
    pub fn to_json(&self) -> Result<String, DatabaseError> {
        Ok(serde_json::to_string_pretty(self)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_database_is_non_empty() {
        let db = MaterialsDatabase::load_builtin();
        assert!(!db.is_empty());
    }

    #[test]
    fn search_by_name_returns_exact_match() {
        let db = MaterialsDatabase::load_builtin();
        let hits = db.search_by_name("AISI 304 (stainless steel)");
        assert!(!hits.is_empty());
    }

    #[test]
    fn search_by_property_returns_subset() {
        let db = MaterialsDatabase::load_builtin();
        let q = PropertyQuery::new(Property::YoungsModulusGPa)
            .min(150.0)
            .max(300.0);
        let hits = db.search_by_property(q);
        for h in &hits {
            assert!(h.mechanical.youngs_modulus_gpa >= 150.0);
            assert!(h.mechanical.youngs_modulus_gpa <= 300.0);
        }
        assert!(!hits.is_empty());
    }

    #[test]
    fn json_round_trip() {
        let db = MaterialsDatabase::load_builtin();
        let s = db.to_json().unwrap();
        let db2 = MaterialsDatabase::from_json(&s).unwrap();
        assert_eq!(db.materials.len(), db2.materials.len());
    }
}