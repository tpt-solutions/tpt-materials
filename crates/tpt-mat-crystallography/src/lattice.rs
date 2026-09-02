//! Lattice parameters for a crystal structure.

use serde::{Deserialize, Serialize};

use tpt_math_linalg_fixed::Vec3;

/// Lattice parameters for a (possibly non-cubic) unit cell.
///
/// Always six components: the three edge lengths `(a, b, c)` and the
/// three angles `(α, β, γ)` in radians.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LatticeParameters {
    /// Lattice parameter `a` (m).
    pub a: f64,
    /// Lattice parameter `b` (m).
    pub b: f64,
    /// Lattice parameter `c` (m).
    pub c: f64,
    /// Angle `α` between `b` and `c` (rad).
    pub alpha: f64,
    /// Angle `β` between `a` and `c` (rad).
    pub beta: f64,
    /// Angle `γ` between `a` and `b` (rad).
    pub gamma: f64,
}

impl LatticeParameters {
    /// Construct a cubic lattice with edge `a` (Å).
    pub fn cubic(a: f64) -> Self {
        Self {
            a,
            b: a,
            c: a,
            alpha: std::f64::consts::FRAC_PI_2,
            beta: std::f64::consts::FRAC_PI_2,
            gamma: std::f64::consts::FRAC_PI_2,
        }
    }

    /// Construct an HCP lattice with `a` and `c` parameters (m).
    pub fn hcp(a: f64, c: f64) -> Self {
        Self {
            a,
            b: a,
            c,
            alpha: std::f64::consts::FRAC_PI_2,
            beta: std::f64::consts::FRAC_PI_2,
            gamma: std::f64::consts::FRAC_PI_2 * 2.0 / 3.0,
        }
    }

    /// Lattice vectors in the sample (cartesian) frame for the canonical
    /// cell orientation (`a` along x, `b` in the x-y plane).
    ///
    /// The matrix is `[a₁ a₂ a₃]` (column-major in `data`).
    pub fn lattice_vectors(&self) -> [Vec3; 3] {
        let a = self.a;
        let b = self.b;
        let c = self.c;
        let ca = self.alpha.cos();
        let cb = self.beta.cos();
        let cg = self.gamma.cos();
        let sg = self.gamma.sin();
        let a1 = Vec3::new(a, 0.0, 0.0);
        let a2 = Vec3::new(b * cg, b * sg, 0.0);
        let a3 = Vec3::new(
            c * cb,
            c * (ca - cb * cg) / sg,
            c * (1.0 - cb * cb - ((ca - cb * cg) / sg).powi(2)).sqrt(),
        );
        [a1, a2, a3]
    }
}
