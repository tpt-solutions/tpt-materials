//! First-order Voigt / Reuss averages for the effective elastic stiffness
//! of an `N`-phase composite.

use tpt_math_linalg_fixed::Vec6;

use tpt_mat_crystal_plasticity::SymmetricFourthOrder;

/// Voigt (iso-strain) average of `N` stiffness tensors:
///
/// `C_V = Σ_r f_r C_r`
///
/// Strictly an upper bound for `C_eff`.  Returns the zero stiffness
/// for an empty list.  `f_r` need not sum to 1; the result is
/// implicitly normalised.
pub fn voigt(phases: &[(SymmetricFourthOrder, f64)]) -> SymmetricFourthOrder {
    let total_f: f64 = phases.iter().map(|(_, f)| *f).sum();
    if total_f <= 0.0 || phases.is_empty() {
        return SymmetricFourthOrder::new([[0.0_f64; 6]; 6]);
    }
    let mut acc = [[0.0_f64; 6]; 6];
    for (c, f) in phases {
        let w = f / total_f;
        for i in 0..6 {
            for j in 0..6 {
                acc[i][j] += w * c.data[i][j];
            }
        }
    }
    SymmetricFourthOrder::new(acc)
}

/// Reuss (iso-stress) average of `N` compliances:
///
/// `S_R = Σ_r f_r S_r = Σ_r f_r C_r^{-1}`
///
/// Strictly a lower bound for the effective stiffness; returned as a
/// 6x6 compliance which the caller can invert if a stiffness is
/// required.
pub fn reuss(phases: &[(SymmetricFourthOrder, f64)]) -> [[f64; 6]; 6] {
    let total_f: f64 = phases.iter().map(|(_, f)| *f).sum();
    if total_f <= 0.0 || phases.is_empty() {
        return [[0.0_f64; 6]; 6];
    }
    let mut acc = [[0.0_f64; 6]; 6];
    for (c, f) in phases {
        let w = f / total_f;
        let s = c.compliance();
        for i in 0..6 {
            for j in 0..6 {
                acc[i][j] += w * s[i][j];
            }
        }
    }
    acc
}

/// Hill (Voigt–Reuss–Hill) arithmetic average of `M_V` and `M_R`.
///
/// `M_VRH = ½ (M_V + M_R)`
///
/// Often used as a quick polycrystal-modulus estimate.  `m_voigt` and
/// `m_reuss` are scalar moduli (e.g. Young's modulus, bulk modulus,
/// shear modulus) in the same units.
pub fn voigt_reuss_average(m_voigt: f64, m_reuss: f64) -> f64 {
    0.5 * (m_voigt + m_reuss)
}

/// Voigt/Reuss bounds on the effective Young's modulus and bulk
/// modulus for an isotropic aggregate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VoigtReussBounds {
    /// Bulk modulus upper bound (Voigt).
    pub k_voigt: f64,
    /// Bulk modulus lower bound (Reuss).
    pub k_reuss: f64,
    /// Shear modulus upper bound (Voigt).
    pub g_voigt: f64,
    /// Shear modulus lower bound (Reuss).
    pub g_reuss: f64,
}

impl VoigtReussBounds {
    /// Hill average bulk modulus.
    pub fn k_vrh(&self) -> f64 {
        voigt_reuss_average(self.k_voigt, self.k_reuss)
    }
    /// Hill average shear modulus.
    pub fn g_vrh(&self) -> f64 {
        voigt_reuss_average(self.g_voigt, self.g_reuss)
    }
    /// Convert `K_VRH, G_VRH` to Young's modulus via the isotropic
    /// relation `E = 9 K G / (3 K + G)`.
    pub fn e_vrh(&self) -> f64 {
        let k = self.k_vrh();
        let g = self.g_vrh();
        9.0 * k * g / (3.0 * k + g)
    }
}

/// Compute the Voigt / Reuss bulk and shear bounds for an aggregate of
/// isotropic phases (each phase carries its own `K_r, G_r`).
pub fn voigt_reuss_bounds(phases: &[(f64, f64, f64)]) -> VoigtReussBounds {
    let total_f: f64 = phases.iter().map(|(_, _, f)| *f).sum();
    if total_f <= 0.0 || phases.is_empty() {
        return VoigtReussBounds {
            k_voigt: 0.0,
            k_reuss: 0.0,
            g_voigt: 0.0,
            g_reuss: 0.0,
        };
    }
    let mut k_v = 0.0_f64;
    let mut g_v = 0.0_f64;
    let mut inv_k_r = 0.0_f64;
    let mut inv_g_r = 0.0_f64;
    for &(k, g, f) in phases {
        let w = f / total_f;
        k_v += w * k;
        g_v += w * g;
        inv_k_r += w / k;
        inv_g_r += w / g;
    }
    VoigtReussBounds {
        k_voigt: k_v,
        k_reuss: if inv_k_r > 0.0 { 1.0 / inv_k_r } else { 0.0 },
        g_voigt: g_v,
        g_reuss: if inv_g_r > 0.0 { 1.0 / inv_g_r } else { 0.0 },
    }
}

/// `E` from isotropic `K, G` via `E = 9 K G / (3 K + G)`.
pub fn youngs_from_k_g(k: f64, g: f64) -> f64 {
    9.0 * k * g / (3.0 * k + g)
}

/// Poisson's ratio `ν` from isotropic `K, G` via `ν = (3 K - 2 G) / (2 (3 K + G))`.
pub fn poisson_from_k_g(k: f64, g: f64) -> f64 {
    (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g))
}

/// Helper: bulk modulus `K` from isotropic `E, ν` via `K = E / (3 (1 - 2 ν))`.
pub fn k_from_e_nu(e: f64, nu: f64) -> f64 {
    e / (3.0 * (1.0 - 2.0 * nu))
}

/// Helper: shear modulus `G` from isotropic `E, ν` via `G = E / (2 (1 + ν))`.
pub fn g_from_e_nu(e: f64, nu: f64) -> f64 {
    e / (2.0 * (1.0 + nu))
}

/// Helper: uniaxial `σ_xx` for strain `E_xx` with isotropic `E, ν`
/// and zero lateral strain.  Currently unused outside tests but
/// retained as a sanity helper.
pub fn _uniaxial_stress_xx(e: f64, nu: f64, eps_xx: f64) -> Vec6 {
    let c = SymmetricFourthOrder::isotropic(e, nu);
    c.contract(Vec6::new(eps_xx, 0.0, 0.0, 0.0, 0.0, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn k_g_from_e_nu(e: f64, nu: f64) -> (f64, f64) {
        (k_from_e_nu(e, nu), g_from_e_nu(e, nu))
    }

    #[test]
    fn voigt_is_upper_bound_for_two_phase_mix() {
        // Soft + stiff phases: Voigt >= Reuss.
        let (k_s, g_s) = k_g_from_e_nu(70_000.0, 0.33); // Al-like
        let (k_h, g_h) = k_g_from_e_nu(400_000.0, 0.25); // Steel-like
        let phases = [(k_s, g_s, 0.5), (k_h, g_h, 0.5)];
        let b = voigt_reuss_bounds(&phases);
        assert!(b.k_voigt >= b.k_reuss);
        assert!(b.g_voigt >= b.g_reuss);
    }

    #[test]
    fn voigt_reuss_single_phase_reproduces_input() {
        let (k, g) = k_g_from_e_nu(200_000.0, 0.3);
        let phases = [(k, g, 1.0)];
        let b = voigt_reuss_bounds(&phases);
        assert!((b.k_voigt - k).abs() < 1e-6);
        assert!((b.k_reuss - k).abs() < 1e-6);
        assert!((b.g_voigt - g).abs() < 1e-6);
        assert!((b.g_reuss - g).abs() < 1e-6);
    }

    #[test]
    fn e_vrh_recovers_youngs_modulus_for_single_phase() {
        let e = 200_000.0_f64;
        let nu = 0.3_f64;
        let (k, g) = k_g_from_e_nu(e, nu);
        let phases = [(k, g, 1.0)];
        let b = voigt_reuss_bounds(&phases);
        let e_calc = b.e_vrh();
        assert!((e_calc - e).abs() < 1.0, "got {e_calc}, expected {e}");
    }

    #[test]
    fn voigt_stiffness_average_matches_6x6() {
        let c_a = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c_b = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let c = voigt(&[(c_a.clone(), 0.3), (c_b.clone(), 0.7)]);
        // C_xx,xx should be 0.3 * c_a[0][0] + 0.7 * c_b[0][0]
        let expected = 0.3 * c_a.data[0][0] + 0.7 * c_b.data[0][0];
        assert!((c.data[0][0] - expected).abs() < 1e-6);
    }

    #[test]
    fn poisson_from_k_g_returns_nu() {
        let e = 200_000.0;
        let nu = 0.3;
        let (k, g) = k_g_from_e_nu(e, nu);
        let nu2 = poisson_from_k_g(k, g);
        assert!((nu2 - nu).abs() < 1e-6);
    }
}
