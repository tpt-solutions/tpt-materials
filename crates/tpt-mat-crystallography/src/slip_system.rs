//! Slip system, Schmid tensor, and resolved shear stress.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use tpt_math_linalg_fixed::{SymMat3, Vec3, Vec6};

use crate::MillerIndex;

/// Family label for a slip system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SlipFamily {
    /// FCC `{111}⟨110⟩`.
    Fcc110,
    /// BCC `{110}⟨111⟩`.
    Bcc111,
    /// HCP basal `(0001)⟨11-20⟩`.
    HcpBasal,
    /// HCP prismatic `{10-10}⟨11-20⟩`.
    HcpPrismatic,
    /// HCP pyramidal ⟨a⟩ `{10-11}⟨11-20⟩`.
    HcpPyramidalA,
    /// HCP pyramidal ⟨c+a⟩ `{11-22}⟨11-23⟩`.
    HcpPyramidalCA,
    /// User-supplied family.
    Other,
}

/// A single slip system.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SlipSystem {
    /// Miller index of the slip plane (e.g. `(1 1 1)` for FCC octahedral).
    pub plane: MillerIndex,
    /// Miller index of the slip direction in 3-index form.
    pub slip_direction_miller: [i32; 3],
    /// Unit-normalised plane normal in the crystal frame.
    pub plane_normal: Vec3,
    /// Unit-normalised slip direction in the crystal frame.
    pub slip_direction: Vec3,
    /// Critical resolved shear stress (CRSS) at reference conditions.
    /// Units: MPa.  Phase 2 wires this to the hardening law.
    pub critical_resolved_shear_stress: f64,
    /// Family label.
    pub family: SlipFamily,
}

/// Errors returned by [`SlipSystem`] operations.
#[derive(Debug, Error, PartialEq)]
pub enum SlipSystemError {
    /// Slip direction was zero.
    #[error("slip direction must be a non-zero vector")]
    ZeroDirection,
    /// Plane normal was zero.
    #[error("plane normal must be a non-zero vector")]
    ZeroNormal,
}

impl SlipSystem {
    /// Construct a slip system from arbitrary non-zero `direction` and
    /// `normal`.  Both are normalised internally.
    pub fn from_vectors(
        direction: Vec3,
        normal: Vec3,
        critical_resolved_shear_stress: f64,
        family: SlipFamily,
    ) -> Result<Self, SlipSystemError> {
        if direction.norm_sq() == 0.0 {
            return Err(SlipSystemError::ZeroDirection);
        }
        if normal.norm_sq() == 0.0 {
            return Err(SlipSystemError::ZeroNormal);
        }
        Ok(Self {
            plane: MillerIndex::new(0, 0, 0),
            slip_direction_miller: [0, 0, 0],
            plane_normal: normal.normalized(),
            slip_direction: direction.normalized(),
            critical_resolved_shear_stress,
            family,
        })
    }

    /// Schmid tensor `P = sym(s ⊗ n) = ½(s ⊗ n + n ⊗ s)`.
    ///
    /// For a stress `σ`, `σ : P` gives the resolved shear stress on the
    /// slip system.
    pub fn schmid_tensor(slip_direction: Vec3, plane_normal: Vec3) -> SymMat3 {
        slip_direction.sym_outer(plane_normal)
    }

    /// Resolved shear stress `τ = σ : P` for stress `σ` in Voigt form.
    ///
    /// `sigma` is a [`Vec6`] in the same frame as `slip_direction` and
    /// `plane_normal`.  Engineering-shear convention is used (matches
    /// `SymMat3::double_dot`).
    pub fn resolved_shear_stress(sigma: Vec6, slip_direction: Vec3, plane_normal: Vec3) -> f64 {
        let p = Self::schmid_tensor(slip_direction, plane_normal);
        sigma.double_dot(Vec6::from_sym_mat3(p))
    }

    /// Max Schmid factor for a tensile axis `t` across a list of slip
    /// systems.  Assumes the slip systems are already in the sample
    /// frame (or, equivalently, `t` is rotated into the crystal
    /// frame).  Returns 0 for an empty list.
    pub fn max_schmid_factor(slips: &[SlipSystem], tensile_axis: [f64; 3]) -> f64 {
        if slips.is_empty() {
            return 0.0;
        }
        let t = Vec3::from(tensile_axis).normalized();
        slips
            .iter()
            .map(|s| {
                let m = (s.slip_direction.dot(t)) * (s.plane_normal.dot(t));
                m.abs()
            })
            .fold(0.0_f64, f64::max)
    }
}

// (helper trait removed — use `Vec6::from_sym_mat3` directly)

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn schmid_tensor_is_symmetric() {
        let s = Vec3::new(1.0, 0.0, 0.0);
        let n = Vec3::new(0.0, 0.0, 1.0);
        let p = SlipSystem::schmid_tensor(s, n);
        // off-diagonal symmetry:
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(p.get(i, j), p.get(j, i), epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn rss_matches_definition_for_simple_stress() {
        // σ = τ (e_x ⊗ e_z + e_z ⊗ e_x); then τ on (s=e_x, n=e_z) == σ_xz.
        let s = Vec3::new(1.0, 0.0, 0.0);
        let n = Vec3::new(0.0, 0.0, 1.0);
        let sigma = Vec6::new(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let tau = SlipSystem::resolved_shear_stress(sigma, s, n);
        assert_relative_eq!(tau, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn max_schmid_factor_matches_theory_for_fcc() {
        // FCC under tension along [001]: max |cos φ cos λ| ≈ 0.4082
        // (Schmid factor of (1 1 -1)[1 1 0] under e_z tension).
        let slips = crate::CrystalStructure::FCC.slip_systems();
        let m = SlipSystem::max_schmid_factor(&slips, [0.0, 0.0, 1.0]);
        assert!((m - 0.4082).abs() < 1e-3, "got {m}");
    }
}
