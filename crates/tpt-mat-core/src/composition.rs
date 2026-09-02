//! Chemical composition of a [`Phase`](crate::Phase).
//!
//! A composition is a map from element symbol to fraction.  The basis —
//! atomic, weight, or mole fraction — is tracked explicitly because the
//! three are not interchangeable.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Basis in which `Composition::fractions` are stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CompositionBasis {
    /// Atomic fraction (count of atoms / total atoms).  Sum = 1.
    Atomic,
    /// Weight (mass) fraction.  Sum = 1.
    Weight,
    /// Mole fraction.  Sum = 1.
    Mole,
}

/// Errors returned by [`Composition`] operations.
#[derive(Debug, Error, PartialEq)]
pub enum CompositionError {
    /// Fractions did not sum to within tolerance of 1.
    #[error("fractions sum to {0}, expected 1.0 within tolerance {1}")]
    NotUnit(f64, f64),
    /// Composition was empty.
    #[error("composition must contain at least one element")]
    Empty,
    /// An element symbol was not found in the periodic table.
    #[error("unknown element symbol: {0:?}")]
    UnknownElement(String),
}

/// Chemical composition of a phase.
///
/// Fractions are non-negative and (within `tol`) sum to 1 in the chosen
/// basis.  Construction normalises and validates.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Composition {
    basis: CompositionBasis,
    /// Element symbol → fraction.  Sorted by symbol via `BTreeMap` for
    /// deterministic serialisation.
    fractions: BTreeMap<String, f64>,
}

impl Composition {
    /// Build a composition, normalising to unit sum.
    ///
    /// `tol` is the absolute deviation from 1 that is accepted in the
    /// _post-normalisation_ sum.  Empty maps are rejected.
    pub fn new(
        basis: CompositionBasis,
        mut fractions: BTreeMap<String, f64>,
        tol: f64,
    ) -> Result<Self, CompositionError> {
        if fractions.is_empty() {
            return Err(CompositionError::Empty);
        }
        let sum: f64 = fractions.values().copied().sum();
        if sum <= 0.0 {
            return Err(CompositionError::NotUnit(sum, tol));
        }
        for v in fractions.values_mut() {
            *v /= sum;
        }
        let new_sum: f64 = fractions.values().sum();
        if (new_sum - 1.0).abs() > tol {
            return Err(CompositionError::NotUnit(new_sum, tol));
        }
        Ok(Self { basis, fractions })
    }

    /// Basis of the stored fractions.
    pub fn basis(&self) -> CompositionBasis {
        self.basis
    }

    /// Fraction of `element` in the stored basis.  Returns `None` if the
    /// element is not present.
    pub fn fraction(&self, element: &str) -> Option<f64> {
        self.fractions.get(element).copied()
    }

    /// All `(element, fraction)` pairs in lexicographic order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, f64)> {
        self.fractions.iter().map(|(k, v)| (k.as_str(), *v))
    }

    /// Convert to a different basis using atomic masses from the
    /// built-in periodic table.
    ///
    /// The conversion goes through mole fractions internally so it is
    /// exact (within numerical tolerance) for both Atomic↔Mole and
    /// Weight↔Mole directions.
    pub fn convert(&self, target: CompositionBasis) -> Result<Self, CompositionError> {
        if target == self.basis {
            return Ok(self.clone());
        }
        // Convert current -> Mole (always possible if elements are known).
        let moles = self.to_moles()?;
        Ok(match target {
            CompositionBasis::Mole | CompositionBasis::Atomic => moles,
            CompositionBasis::Weight => moles.to_weight()?,
        })
    }

    fn to_moles(&self) -> Result<Self, CompositionError> {
        if matches!(self.basis, CompositionBasis::Mole) {
            return Ok(self.clone());
        }
        if matches!(self.basis, CompositionBasis::Atomic) {
            // Atomic fraction ≡ mole fraction for single-phase multi-element.
            return Ok(self.clone());
        }
        // Weight -> Mole: x_i = w_i / M_i / sum(w_j / M_j).
        let mut moles = BTreeMap::new();
        let mut sum = 0.0;
        for (el, w) in self.fractions.iter() {
            let m = tpt_mat_constants::PeriodicTable::atomic_mass(el)
                .ok_or_else(|| CompositionError::UnknownElement(el.clone()))?;
            let mi = w / m;
            moles.insert(el.clone(), mi);
            sum += mi;
        }
        for v in moles.values_mut() {
            *v /= sum;
        }
        Ok(Self {
            basis: CompositionBasis::Mole,
            fractions: moles,
        })
    }

    fn to_weight(&self) -> Result<Self, CompositionError> {
        let mut weights = BTreeMap::new();
        let mut sum = 0.0;
        for (el, x) in self.fractions.iter() {
            let m = tpt_mat_constants::PeriodicTable::atomic_mass(el)
                .ok_or_else(|| CompositionError::UnknownElement(el.clone()))?;
            let wi = x * m;
            weights.insert(el.clone(), wi);
            sum += wi;
        }
        for v in weights.values_mut() {
            *v /= sum;
        }
        Ok(Self {
            basis: CompositionBasis::Weight,
            fractions: weights,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fe_cr() -> BTreeMap<String, f64> {
        let mut m = BTreeMap::new();
        m.insert("Fe".into(), 0.7);
        m.insert("Cr".into(), 0.3);
        m
    }

    #[test]
    fn construct_normalises() {
        let c = Composition::new(CompositionBasis::Atomic, fe_cr(), 1e-12).unwrap();
        let sum: f64 = c.iter().map(|(_, v)| v).sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn empty_is_rejected() {
        assert!(matches!(
            Composition::new(CompositionBasis::Atomic, BTreeMap::new(), 1e-12),
            Err(CompositionError::Empty)
        ));
    }

    #[test]
    fn weight_to_mole_roundtrip() {
        let w = Composition::new(CompositionBasis::Weight, fe_cr(), 1e-12).unwrap();
        let m = w.convert(CompositionBasis::Mole).unwrap();
        let back = m.convert(CompositionBasis::Weight).unwrap();
        for (el, f) in w.iter() {
            assert!((f - back.fraction(el).unwrap()).abs() < 1e-10);
        }
    }
}
