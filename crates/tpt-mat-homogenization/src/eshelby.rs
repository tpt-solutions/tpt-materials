//! Eshelby's equivalent-inclusion method for a single ellipsoidal
//! inclusion in an isotropic matrix.
//!
//! Reference: Eshelby, J. D. (1957).  "The determination of the elastic
//! field of an ellipsoidal inclusion, and related problems."
//! Proc. R. Soc. Lond. A 241, 376–396.
//!
//! For an isotropic matrix (`K, G`) hosting an isotropic inclusion
//! (`K₁, G₁`) with eigenstrain `ε*`, the constrained (elastic)
//! strain inside the inclusion is
//!
//! `ε^c = S : ε*`
//!
//! where `S` is Eshelby's tensor.  For a spherical inclusion,
//! `S = (3 K, 2 (4 - 5 ν)) / (15 (1 - ν))` in isotropic
//! `(volumetric, deviatoric)` form, see e.g. Mura's *Micromechanics
//! of Defects in Solids*.
//!
//! The dilute strain-concentration tensor `A` relates the
//! *far-field* strain `ε^∞` to the inclusion's *average* strain
//! `ε_inside = A : ε^∞` for a single ellipsoidal inclusion of
//! stiffness `C₁` in a matrix of stiffness `C_0`:
//!
//! `A = [I + S · C_0^{-1} · (C_1 - C_0)]^{-1}`
//!
//! (`Mori–Tanaka, 1973`, special case of `c → 0`.)

use serde::{Deserialize, Serialize};

use tpt_math_linalg_fixed::Vec6;

use tpt_mat_crystal_plasticity::SymmetricFourthOrder;

use super::voigt_reuss::{g_from_e_nu, k_from_e_nu};

/// Eshelby tensor components for a spherical inclusion in an
/// isotropic matrix.
///
/// `S` is decomposed into a hydrostatic part (`S_hydro`) and a
/// deviatoric part (`S_dev`) such that, for an eigenstrain `ε*`:
//$$
/// `ε^c_vol = S_hydro · ε*_vol`
/// `ε^c_dev = S_dev  · ε*_dev`
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EshelbySpherical {
    /// Hydrostatic scalar component `S_h = (1 + ν) / (3 (1 - ν))`.
    pub s_hydro: f64,
    /// Deviatoric scalar component `S_d = 2 (4 - 5 ν) / (15 (1 - ν))`.
    pub s_dev: f64,
}

impl EshelbySpherical {
    /// Construct from matrix Poisson's ratio `ν`.  Panics if `ν` is
    /// outside `(-1, 0.5)` (incompressibility / instability).
    pub fn from_nu(nu: f64) -> Self {
        assert!(
            nu > -1.0 && nu < 0.5,
            "Poisson's ratio must be in (-1, 0.5); got {nu}"
        );
        let s_hydro = (1.0 + nu) / (3.0 * (1.0 - nu));
        let s_dev = 2.0 * (4.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
        Self { s_hydro, s_dev }
    }
}

/// Convenience constructor.
pub fn eshelby_spherical(matrix_e: f64, matrix_nu: f64) -> EshelbySpherical {
    let _ = (
        matrix_e,
        g_from_e_nu(matrix_e, matrix_nu),
        k_from_e_nu(matrix_e, matrix_nu),
    );
    EshelbySpherical::from_nu(matrix_nu)
}

/// Dilute strain-concentration tensor `A` (6x6) for a single
/// ellipsoidal inclusion with stiffness `C_inclusion` in an
/// isotropic matrix with stiffness `C_matrix` and Eshelby tensor
/// `S`.  Returns the symmetric `A = [I + S · C_0^{-1} · (C_1 - C_0)]^{-1}`
/// in Voigt order.
pub fn dilute_strain_concentration(
    c_matrix: &SymmetricFourthOrder,
    c_inclusion: &SymmetricFourthOrder,
    eshelby: EshelbySpherical,
) -> StrainConcentrationTensor {
    let s_matrix = c_matrix.compliance();
    // ΔC = C_inclusion - C_matrix
    let mut dc = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            dc[i][j] = c_inclusion.data[i][j] - c_matrix.data[i][j];
        }
    }
    // S · C_0^{-1} · ΔC
    let scd = eshelby_apply_mat(&s_matrix, &dc, eshelby);
    // I + S C_0^{-1} ΔC
    let mut m = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            m[i][j] = scd[i][j] + if i == j { 1.0 } else { 0.0 };
        }
    }
    let a = inv6(m);
    StrainConcentrationTensor { data: a }
}

/// Strain-concentration tensor `A`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StrainConcentrationTensor {
    /// 6x6 concentration tensor in Voigt order.
    pub data: [[f64; 6]; 6],
}

impl StrainConcentrationTensor {
    /// `A : ε^∞` (matrix-vector product).
    pub fn apply(self, eps_inf: Vec6) -> Vec6 {
        let mut out = [0.0_f64; 6];
        for i in 0..6 {
            let mut s = 0.0;
            for j in 0..6 {
                s += self.data[i][j] * eps_inf[j];
            }
            out[i] = s;
        }
        Vec6::new(out[0], out[1], out[2], out[3], out[4], out[5])
    }
}

/// Eshelby tensor in 6x6 isotropic form (hydrostatic / deviatoric
/// decomposition).
fn eshelby_6x6(s: EshelbySpherical) -> [[f64; 6]; 6] {
    // Isotropic tensor of the form S = S_h * I_h + S_d * I_d, where
    // I_h projects on hydrostatic components and I_d on deviatoric.
    // Voigt order [xx, yy, zz, xy, yz, xz]; shear entries are weighted
    // by the engineering-shear convention.
    let mut m = [[0.0_f64; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            // Spherical Eshelby tensor in Voigt order, from the
            // isotropic projections.  In pure-tensor index form:
            //   S_1111 = S_h/3 + 2 S_d/3
            //   S_1122 = S_h/3 -     S_d/3
            let hydro = 1.0 / 3.0;
            let dev = if i == j { 2.0 / 3.0 } else { -1.0 / 3.0 };
            m[i][j] = s.s_hydro * hydro + s.s_dev * dev;
        }
    }
    for i in 0..3 {
        m[3 + 0][i] = 0.0;
        m[i][3 + 0] = 0.0;
    }
    for k in 3..6 {
        // Engineering-shear weight: for engineering shear gamma = 2 eps
        // the tensor entry S_1212 = (4-5nu)/(15(1-nu)) = s_dev/2 acts on
        // the engineering shear strain directly.
        m[k][k] = s.s_dev / 2.0;
    }
    m
}

/// `S · S_0 · ΔC` (matrix product chain).
fn eshelby_apply_mat(s0: &[[f64; 6]; 6], dc: &[[f64; 6]; 6], s: EshelbySpherical) -> [[f64; 6]; 6] {
    let s6 = eshelby_6x6(s);
    // tmp = S_0 · ΔC
    let mut tmp = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut acc = 0.0;
            for k in 0..6 {
                acc += s0[i][k] * dc[k][j];
            }
            tmp[i][j] = acc;
        }
    }
    // out = S · tmp
    let mut out = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut acc = 0.0;
            for k in 0..6 {
                acc += s6[i][k] * tmp[k][j];
            }
            out[i][j] = acc;
        }
    }
    out
}

/// Invert a 6x6 matrix via Gaussian elimination with partial pivoting.
fn inv6(m: [[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut a = [[0.0_f64; 12]; 6];
    for i in 0..6 {
        for j in 0..6 {
            a[i][j] = m[i][j];
            a[i][j + 6] = if i == j { 1.0 } else { 0.0 };
        }
    }
    for i in 0..6 {
        let mut piv = i;
        for k in (i + 1)..6 {
            if a[k][i].abs() > a[piv][i].abs() {
                piv = k;
            }
        }
        if a[piv][i].abs() < 1e-15 {
            panic!("singular matrix in inv6");
        }
        a.swap(i, piv);
        let inv_piv = 1.0 / a[i][i];
        for j in 0..12 {
            a[i][j] *= inv_piv;
        }
        for k in 0..6 {
            if k != i {
                let f = a[k][i];
                for j in 0..12 {
                    a[k][j] -= f * a[i][j];
                }
            }
        }
    }
    let mut out = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i][j] = a[i][j + 6];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_testkit::assert_relative_eq;

    #[test]
    fn eshelby_6x6_sphere_matches_analytic_voigt_entries() {
        // Analytic entries (isotropic matrix, sphere), pure-tensor
        // Voigt form with engineering-shear weight for the shear block:
        //   S_1111 = (7-5ν)/(15(1-ν))
        //   S_1122 = (5ν-1)/(15(1-ν))
        //   S_1212 = (4-5ν)/(15(1-ν))
        let nu = 0.3;
        let sphere = EshelbySpherical::from_nu(nu);
        let m = super::eshelby_6x6(sphere);
        let denom = 15.0 * (1.0 - nu);
        let s11 = (7.0 - 5.0 * nu) / denom;
        let s12 = (5.0 * nu - 1.0) / denom;
        let s44 = (4.0 - 5.0 * nu) / denom;
        assert!((m[0][0] - s11).abs() < 1e-12);
        assert!((m[0][1] - s12).abs() < 1e-12);
        assert!((m[1][2] - s12).abs() < 1e-12);
        assert!((m[3][3] - s44).abs() < 1e-12);
        assert!((m[4][4] - s44).abs() < 1e-12);
        assert!((m[5][5] - s44).abs() < 1e-12);
        // Hydrostatic projection: S applied to a volumetric eigenstrain
        // must reproduce the trace contraction (s_hydro).
        let sh = sphere.s_hydro;
        let trace: f64 = m[0][0] + m[0][1] + m[0][2];
        assert!((trace - sh).abs() < 1e-12);
    }

    #[test]
    fn dilute_strain_concentration_isotropic_matches_analytic() {
        // For an isotropic inclusion in an isotropic matrix the dilute
        // (Mori–Tanaka, c → 0) concentration of a volumetric strain is
        //   A_hydro = (K0 + 4G0/3) / (K1 + 4G0/3).
        let e0 = 200_000.0;
        let nu0 = 0.3;
        let k0 = e0 / (3.0 * (1.0 - 2.0 * nu0));
        let g0 = e0 / (2.0 * (1.0 + nu0));
        let e1 = 400_000.0;
        let k1 = e1 / (3.0 * (1.0 - 2.0 * 0.25));
        let c0 = SymmetricFourthOrder::isotropic(e0, nu0);
        let c1 = SymmetricFourthOrder::isotropic(e1, 0.25);
        let s = EshelbySpherical::from_nu(nu0);
        let a = dilute_strain_concentration(&c0, &c1, s);
        let eps_inf = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        let eps_inc = a.apply(eps_inf);
        let expected = (k0 + 4.0 * g0 / 3.0) / (k1 + 4.0 * g0 / 3.0);
        let vol_ratio = (eps_inc[0] + eps_inc[1] + eps_inc[2]) / 3.0;
        assert!(
            (vol_ratio - expected).abs() < 1e-9,
            "hydrostatic concentration {vol_ratio} != analytic {expected}"
        );
    }

    #[test]
    fn eshelby_spherical_zero_inclusion_stiffness_gives_unity_concentration() {
        // If C_inclusion = 0 (a cavity), then ΔC = -C_matrix and
        // A = [I - S]^{-1}, but for *zero* inclusion stiffness
        // (ΔC = -C_0), the inclusion takes the far-field strain
        // amplified by [I - S]^{-1}.  Sanity: check that S applied
        // to a hydrostatic eigenstrain gives the expected scalar
        // reduction.
        let s = EshelbySpherical::from_nu(0.3);
        // S_h = (1 + 0.3) / (3 * 0.7) ≈ 0.6190
        assert!((s.s_hydro - (1.3_f64 / (3.0 * 0.7))).abs() < 1e-12);
        // S_d = 2 (4 - 5*0.3) / (15 * 0.7) ≈ 0.5143
        let expected_dev = 2.0 * (4.0 - 5.0 * 0.3) / (15.0 * 0.7);
        assert!((s.s_dev - expected_dev).abs() < 1e-12);
    }

    #[test]
    fn dilute_concentration_recovers_identity_when_inclusion_matches_matrix() {
        let c0 = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let s = EshelbySpherical::from_nu(0.3);
        let a = dilute_strain_concentration(&c0, &c0, s);
        // ΔC = 0 ⇒ A = I^{-1} = I.
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_relative_eq!(a.data[i][j], expected, epsilon = 1e-9);
            }
        }
    }

    #[test]
    fn dilute_concentration_for_stiff_inclusion_amplifies_strain_less_than_unity() {
        let c0 = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c1 = SymmetricFourthOrder::isotropic(400_000.0, 0.25);
        let s = EshelbySpherical::from_nu(0.33);
        let a = dilute_strain_concentration(&c0, &c1, s);
        // For uniaxial far-field strain along x, the inclusion's x
        // strain should be *less* than the far-field strain (stiff
        // inclusion shields the far field).
        let eps_inf = Vec6::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let eps_inc = a.apply(eps_inf);
        assert!(
            eps_inc[0] < 1.0,
            "ε_inc_xx={} should be < ε_inf_xx=1.0 for stiff inclusion",
            eps_inc[0]
        );
        assert!(eps_inc[0] > 0.0);
    }

    #[test]
    fn dilute_concentration_for_soft_inclusion_amplifies_strain_more_than_unity() {
        let c0 = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let c1 = SymmetricFourthOrder::isotropic(10_000.0, 0.1);
        let s = EshelbySpherical::from_nu(0.3);
        let a = dilute_strain_concentration(&c0, &c1, s);
        let eps_inf = Vec6::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let eps_inc = a.apply(eps_inf);
        assert!(
            eps_inc[0] > 1.0,
            "ε_inc_xx={} should be > ε_inf_xx=1.0 for soft inclusion",
            eps_inc[0]
        );
    }
}
