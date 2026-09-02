//! Crystal-plasticity model configuration.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use tpt_mat_crystallography::{CrystalStructure, SlipSystem};
use tpt_mat_hardening::Hardening;

use crate::elastic::SymmetricFourthOrder;

/// Rate-sensitivity parameter for the viscoplastic power-law flow rule.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RateSensitivity {
    /// Reference slip rate `γ̇_0` (1/s).
    pub reference_strain_rate: f64,
    /// Rate-sensitivity exponent `n` (dimensionless).
    pub exponent: f64,
}

impl Default for RateSensitivity {
    fn default() -> Self {
        Self {
            reference_strain_rate: 1.0e-3,
            exponent: 20.0,
        }
    }
}

/// Errors raised when constructing a [`CrystalPlasticityModel`].
#[derive(Debug, Error, PartialEq)]
pub enum CrystalPlasticityModelError {
    /// `elastic_tensor` matrix is not 6x6.
    #[error("elastic tensor must be 6x6")]
    BadElasticSize,
    /// `n_slip` from `crystal_structure.slip_systems()` and the supplied
    /// `slip_systems` vector disagree.
    #[error("slip-systems count {0} does not match crystal structure default {1}")]
    SlipCountMismatch(usize, usize),
}

/// Single-crystal constitutive model.
///
/// Captures the crystal structure, the explicit slip systems, the
/// hardening law, the rate-sensitivity parameters, and the elastic
/// stiffness tensor.  The model is evaluated per integration point
/// inside an FEM element (see [`CpFemSolver`](crate::CpFemSolver)) or in
/// a stand-alone single-point update.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CrystalPlasticityModel {
    /// Crystal structure (FCC, BCC, ...).
    pub crystal_structure: CrystalStructure,
    /// Explicit slip systems (overrides `crystal_structure.slip_systems()`).
    pub slip_systems: Vec<SlipSystem>,
    /// Hardening law.
    pub hardening_law: Hardening,
    /// Viscoplastic flow-rule parameters.
    pub rate_sensitivity: RateSensitivity,
    /// Elastic stiffness in Voigt 6x6 form (MPa).
    pub elastic_tensor: SymmetricFourthOrder,
}

impl CrystalPlasticityModel {
    /// Build from a [`CrystalStructure`]; the slip systems are filled in
    /// from the structure's defaults.
    pub fn from_crystal_structure(
        crystal_structure: CrystalStructure,
        hardening_law: Hardening,
        rate_sensitivity: RateSensitivity,
        elastic_tensor: SymmetricFourthOrder,
    ) -> Result<Self, CrystalPlasticityModelError> {
        let slip_systems = crystal_structure.slip_systems();
        Ok(Self {
            crystal_structure,
            slip_systems,
            hardening_law,
            rate_sensitivity,
            elastic_tensor,
        })
    }

    /// Number of slip systems.
    pub fn n_slip(&self) -> usize {
        self.slip_systems.len()
    }
}