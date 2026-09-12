//! Representative Volume Element (RVE) data model and simple
//! homogenizers.
//!
//! This module *does not* perform a full FEM or FFT-based
//! homogenization; for that, a downstream `tpt-fem` integration is
//! required (deferred to a follow-up phase).  Instead it provides:
//!
//! - [`Rve`]: a list of grains with orientations and stiffnesses.
//! - [`HomogenizationScheme`]: the available analytical bounds.
//! - [`SimpleHomogenizer`]: a convenience driver for Voigt / Reuss /
//!   self-consistent homogenization over the grain list.

use serde::{Deserialize, Serialize};

use tpt_mat_crystal_plasticity::SymmetricFourthOrder;
use tpt_math_linalg_fixed::Vec6;

use crate::homogenization::{dilute_strain_concentration, reuss, voigt, EshelbySpherical};

/// A single grain in the RVE.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RveGrain {
    /// Optional label (e.g. grain ID).
    pub label: String,
    /// Volume fraction.  Should sum to 1 across the RVE for proper
    /// homogenization; the helpers normalise internally if not.
    pub volume_fraction: f64,
    /// Orientation matrix `R` mapping the crystal frame to the sample
    /// frame (`v_sample = R · v_crystal`).
    pub orientation: [[f64; 3]; 3],
    /// Anisotropic elastic stiffness in the crystal frame.
    pub stiffness: SymmetricFourthOrder,
}

impl RveGrain {
    /// Convenience constructor.
    pub fn new(
        label: impl Into<String>,
        volume_fraction: f64,
        orientation: [[f64; 3]; 3],
        stiffness: SymmetricFourthOrder,
    ) -> Self {
        Self {
            label: label.into(),
            volume_fraction,
            orientation,
            stiffness,
        }
    }

    /// Rotate `v` from the crystal frame to the sample frame.
    pub fn rotate_to_sample(&self, v: Vec6) -> Vec6 {
        // Apply the 3x3 rotation to the 3x3 sub-block; engineering
        // shear components transform like `σ'_xy = R_x R_y σ_xy` etc.
        let r = self.orientation;
        let m = v.to_sym_mat3_data();
        let mut out = [0.0_f64; 6];
        for i in 0..3 {
            for j in 0..3 {
                let mut acc = 0.0;
                for k in 0..3 {
                    for l in 0..3 {
                        acc += r[i][k] * r[j][l] * m[k][l];
                    }
                }
                out[voigt_idx(i, j)] = acc;
            }
        }
        Vec6::new(out[0], out[1], out[2], out[3], out[4], out[5])
    }

    /// Stiffness tensor rotated from the crystal frame into the
    /// sample frame.
    ///
    /// The stored 6x6 matrix uses engineering shear strain (`σ = C·E`
    /// with `E_j = (2 - δ) ε`), the standard no-double Voigt layout
    /// whose entries *are* the 4th-order tensor components.  Rotating
    /// the 4th-order tensor by `R` and re-compressing gives the
    /// symmetrised basis `Q[i][k] = Σ_{(m,n)∈pair(k)} R_{a,m} R_{b,n}`
    /// (both orderings of each shear pair), so `C' = Q·C·Qᵀ`
    /// reproduces the exact `R ⊗₄ R` rotation to machine precision
    /// and preserves the tensor invariants and positivity by
    /// construction.
    pub fn rotated_stiffness(&self) -> SymmetricFourthOrder {
        let r = self.orientation;
        let q = |i: usize, k: usize| -> f64 {
            let (a, b) = voigt_pair(i);
            match k {
                0..=2 => r[a][k] * r[b][k],
                3 => r[a][0] * r[b][1] + r[a][1] * r[b][0],
                4 => r[a][1] * r[b][2] + r[a][2] * r[b][1],
                5 => r[a][0] * r[b][2] + r[a][2] * r[b][0],
                _ => unreachable!("Voigt index {k}"),
            }
        };
        let s = &self.stiffness.data;
        let mut c = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let mut acc = 0.0;
                for k in 0..6 {
                    for l in 0..6 {
                        acc += q(i, k) * s[k][l] * q(j, l);
                    }
                }
                c[i][j] = acc;
            }
        }
        SymmetricFourthOrder::new(c)
    }
}

#[inline]
fn voigt_idx(i: usize, j: usize) -> usize {
    match (i, j) {
        (0, 0) => 0,
        (1, 1) => 1,
        (2, 2) => 2,
        (0, 1) | (1, 0) => 3,
        (1, 2) | (2, 1) => 4,
        (0, 2) | (2, 0) => 5,
        _ => panic!("Voigt index out of range"),
    }
}

#[inline]
fn voigt_pair(idx: usize) -> (usize, usize) {
    match idx {
        0 => (0, 0),
        1 => (1, 1),
        2 => (2, 2),
        3 => (0, 1),
        4 => (1, 2),
        5 => (0, 2),
        _ => panic!("Voigt index {idx} out of range"),
    }
}

/// The RVE: a collection of grains.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Rve {
    /// All grains.
    pub grains: Vec<RveGrain>,
}

impl Rve {
    /// Construct from a list of grains.
    pub fn new(grains: Vec<RveGrain>) -> Self {
        Self { grains }
    }

    /// Empty RVE.
    pub fn empty() -> Self {
        Self { grains: Vec::new() }
    }

    /// Total volume fraction (typically 1.0).
    pub fn total_volume_fraction(&self) -> f64 {
        self.grains.iter().map(|g| g.volume_fraction).sum()
    }

    /// Number of grains.
    pub fn n_grains(&self) -> usize {
        self.grains.len()
    }

    /// Effective stiffness under the chosen homogenization scheme
    /// (Voigt, Reuss or one-site self-consistent).
    pub fn homogenize(&self, scheme: HomogenizationScheme) -> SymmetricFourthOrder {
        match scheme {
            HomogenizationScheme::Voigt => {
                let phases: Vec<_> = self
                    .grains
                    .iter()
                    .map(|g| (g.rotated_stiffness(), g.volume_fraction))
                    .collect();
                voigt(&phases)
            }
            HomogenizationScheme::Reuss => {
                // Reuss: harmonic average of compliances, returned
                // as a 6x6 compliance; invert to a stiffness.
                let phases: Vec<_> = self
                    .grains
                    .iter()
                    .map(|g| (g.rotated_stiffness(), g.volume_fraction))
                    .collect();
                let s_avg = reuss(&phases);
                inv6(s_avg)
            }
            HomogenizationScheme::SelfConsistent => {
                // One-site self-consistent (Kröner, 1958; Budiansky &
                // Wu, 1962) scheme.
                //
                // Each grain is treated as a spherical Eshelby
                // inclusion in the *unknown* effective medium `C*`; the
                // strain-concentration tensors
                // `A_i = [I + S(C*, ν*) C*^{-1} (C_i - C*)]^{-1}` are
                // built from the Eshelby tensor of the implicit medium
                // and the medium is updated by enforcing the
                // self-consistency condition
                // `C* = ⟨C_i A_i⟩ ⟨A_i⟩^{-1}`.  Iterated to a fixed
                // point this lies strictly between the Voigt and Reuss
                // bounds (it coincides with the Kröner–Budiansky–Wu
                // estimate for isotropic / weakly anisotropic grains).
                let n = self.grains.len();
                if n == 0 {
                    return SymmetricFourthOrder::new([[0.0_f64; 6]; 6]);
                }
                let total_f: f64 = self.grains.iter().map(|g| g.volume_fraction).sum();
                let phases: Vec<_> = self
                    .grains
                    .iter()
                    .map(|g| (g.rotated_stiffness(), g.volume_fraction))
                    .collect();
                let mut c_eff = voigt(&phases);
                const MAX_ITER: usize = 60;
                const TOL: f64 = 1.0e-10;
                // Under-relaxation for strongly anisotropic grains:
                // the undamped Picard map is not a global contraction
                // and can wander behind the Voigt/Reuss envelope.
                const OMEGA: f64 = 0.6;
                for _ in 0..MAX_ITER {
                    let eshelby = EshelbySpherical::from_nu(effective_poisson(&c_eff));
                    let mut weighted = [[0.0_f64; 6]; 6]; // Σ f_i C_i A_i
                    let mut a_avg = [[0.0_f64; 6]; 6]; // Σ f_i A_i
                    for g in &self.grains {
                        let c_i = g.rotated_stiffness();
                        let f_i = if total_f > 0.0 {
                            g.volume_fraction / total_f
                        } else {
                            0.0
                        };
                        let a = dilute_strain_concentration(&c_eff, &c_i, eshelby).data;
                        let c_a = mat6_mul(&c_i.data, &a);
                        for i in 0..6 {
                            for j in 0..6 {
                                weighted[i][j] += f_i * c_a[i][j];
                                a_avg[i][j] += f_i * a[i][j];
                            }
                        }
                    }
                    let c_new = mat6_mul(&weighted, &inv6(a_avg).data);
                    let mut blended = [[0.0_f64; 6]; 6];
                    for i in 0..6 {
                        for j in 0..6 {
                            blended[i][j] =
                                c_eff.data[i][j] + OMEGA * (c_new[i][j] - c_eff.data[i][j]);
                        }
                    }
                    let mut max_diff = 0.0_f64;
                    for i in 0..6 {
                        for j in 0..6 {
                            max_diff = max_diff.max((blended[i][j] - c_eff.data[i][j]).abs());
                        }
                    }
                    c_eff = SymmetricFourthOrder::new(blended);
                    if max_diff < TOL * c_eff.data[0][0].max(1.0) {
                        break;
                    }
                }
                c_eff
            }
        }
    }
}

/// 6x6 matrix product.
fn mat6_mul(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut out = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut acc = 0.0;
            for k in 0..6 {
                acc += a[i][k] * b[k][j];
            }
            out[i][j] = acc;
        }
    }
    out
}

/// Effective-medium Poisson's ratio from the isotropic part of `C`.
fn effective_poisson(c: &SymmetricFourthOrder) -> f64 {
    let d = c.data;
    let c11 = (d[0][0] + d[1][1] + d[2][2]) / 3.0;
    let c12 = (d[0][1] + d[0][2] + d[1][2]) / 3.0;
    let c44 = (d[3][3] + d[4][4] + d[5][5]) / 3.0;
    let k = (c11 + 2.0 * c12) / 3.0;
    let g = c44;
    // ν = (3K − 2G) / (2 (3K + G)).
    let denom = 2.0 * (3.0 * k + g);
    if denom < 1e-15 {
        return 0.3;
    }
    let nu = (3.0 * k - 2.0 * g) / denom;
    nu.clamp(-0.99, 0.499)
}

/// Available homogenization schemes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum HomogenizationScheme {
    /// Voigt (iso-strain) upper bound.
    Voigt,
    /// Reuss (iso-stress) lower bound.
    Reuss,
    /// One-site self-consistent mean field (Kröner / Budiansky–Wu).
    SelfConsistent,
}

/// Aggregate statistics of the RVE.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RveStats {
    /// Number of grains.
    pub n_grains: usize,
    /// Total volume fraction.
    pub total_volume_fraction: f64,
    /// Number of unique orientations (rounded to a coarse grid).
    pub unique_orientations: usize,
}

/// Quick stats.
impl Rve {
    pub fn stats(&self) -> RveStats {
        // Count near-identical orientations (entries within 1°).
        let mut count = 0;
        for (i, a) in self.grains.iter().enumerate() {
            let mut is_unique = true;
            for b in self.grains.iter().take(i) {
                if orientation_close(&a.orientation, &b.orientation, 1.0_f64.to_radians()) {
                    is_unique = false;
                    break;
                }
            }
            if is_unique {
                count += 1;
            }
        }
        RveStats {
            n_grains: self.grains.len(),
            total_volume_fraction: self.total_volume_fraction(),
            unique_orientations: count,
        }
    }
}

fn orientation_close(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3], tol: f64) -> bool {
    let mut max = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            max = max.max((a[i][j] - b[i][j]).abs());
        }
    }
    max < tol
}

/// Quick helper for tests/examples: a simple isotropic homogenizer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct SimpleHomogenizer;

impl SimpleHomogenizer {
    /// Voigt average stiffness (delegates to [`voigt`]).
    pub fn voigt(rve: &Rve) -> SymmetricFourthOrder {
        rve.homogenize(HomogenizationScheme::Voigt)
    }
    /// Reuss average stiffness (delegates to [`reuss`] + invert).
    pub fn reuss(rve: &Rve) -> SymmetricFourthOrder {
        rve.homogenize(HomogenizationScheme::Reuss)
    }
}

fn inv6(m: [[f64; 6]; 6]) -> SymmetricFourthOrder {
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
            panic!("singular");
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
    SymmetricFourthOrder::new(out)
}

// Bring helper trait into scope for `Vec6::to_sym_mat3_data`.
trait SymMat3Ext {
    fn to_sym_mat3_data(self) -> [[f64; 3]; 3];
}
impl SymMat3Ext for Vec6 {
    fn to_sym_mat3_data(self) -> [[f64; 3]; 3] {
        let d = self.data;
        [[d[0], d[3], d[5]], [d[3], d[1], d[4]], [d[5], d[4], d[2]]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity() -> [[f64; 3]; 3] {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    #[test]
    fn single_grain_voigt_recovers_input() {
        let c = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let rve = Rve::new(vec![RveGrain::new("g0", 1.0, identity(), c.clone())]);
        let c_avg = rve.homogenize(HomogenizationScheme::Voigt);
        for i in 0..6 {
            for j in 0..6 {
                assert!((c_avg.data[i][j] - c.data[i][j]).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn empty_rve_yields_zero_stiffness() {
        let rve = Rve::empty();
        let c = rve.homogenize(HomogenizationScheme::Voigt);
        assert_eq!(c.data[0][0], 0.0);
    }

    #[test]
    fn two_grain_voigt_average_is_weighted_mean() {
        let c_a = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c_b = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let rve = Rve::new(vec![
            RveGrain::new("a", 0.3, identity(), c_a.clone()),
            RveGrain::new("b", 0.7, identity(), c_b.clone()),
        ]);
        let c_avg = rve.homogenize(HomogenizationScheme::Voigt);
        let expected = 0.3 * c_a.data[0][0] + 0.7 * c_b.data[0][0];
        assert!((c_avg.data[0][0] - expected).abs() < 1e-6);
    }

    #[test]
    fn self_consistent_lies_strictly_between_reuss_and_voigt() {
        let c_a = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c_b = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let rve = Rve::new(vec![
            RveGrain::new("a", 0.3, identity(), c_a.clone()),
            RveGrain::new("b", 0.7, identity(), c_b.clone()),
        ]);
        let c_v = rve.homogenize(HomogenizationScheme::Voigt);
        let c_r = rve.homogenize(HomogenizationScheme::Reuss);
        let c_sc = rve.homogenize(HomogenizationScheme::SelfConsistent);
        // Bulk and shear components of the SC estimate must lie
        // strictly inside the Voigt/Reuss interval.
        let bulk = |c: &SymmetricFourthOrder| (c.data[0][0] + 2.0 * c.data[0][1]) / 3.0;
        let shear = |c: &SymmetricFourthOrder| c.data[3][3];
        assert!(bulk(&c_r) < bulk(&c_sc) && bulk(&c_sc) < bulk(&c_v));
        assert!(shear(&c_r) < shear(&c_sc) && shear(&c_sc) < shear(&c_v));
        // Symmetry of the returned tensor (new() symmetrises).
        for i in 0..6 {
            for j in 0..6 {
                assert!((c_sc.data[i][j] - c_sc.data[j][i]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn self_consistent_anisotropic_random_grains_stays_between_bounds() {
        // Regression: the undamped Picard map could wander behind the
        // Voigt/Reuss envelope for randomly oriented cubic grains and
        // its fixed point exploded.  With the isotropic Eshelby tensor
        // fixed and an under-relaxed update, the estimate must stay
        // between the bounds and the bulk modulus must match the
        // (orientation-independent) crystal value.
        let mut seed = 7u64;
        let stiffness = SymmetricFourthOrder::cubic(168.4e9, 121.4e9, 75.4e9);
        let grains: Vec<RveGrain> = (0..80)
            .map(|i| {
                RveGrain::new(
                    format!("g{i}"),
                    1.0 / 80.0,
                    random_bunge(&mut seed),
                    stiffness.clone(),
                )
            })
            .collect();
        let rve = Rve::new(grains);
        let c_v = rve.homogenize(HomogenizationScheme::Voigt);
        let c_r = rve.homogenize(HomogenizationScheme::Reuss);
        let c_sc = rve.homogenize(HomogenizationScheme::SelfConsistent);
        let fcc_bulk = (168.4e9 + 2.0 * 121.4e9) / 3.0;
        // Rotation-invariant bulk (valid for fully anisotropic tensors,
        // unlike (C11 + 2 C12) / 3 which only applies to cubes).
        let bulk = |c: &SymmetricFourthOrder| {
            (c.data[0][0]
                + c.data[1][1]
                + c.data[2][2]
                + 2.0 * (c.data[0][1] + c.data[0][2] + c.data[1][2]))
                / 9.0
        };
        let shear = |c: &SymmetricFourthOrder| c.data[3][3];
        assert!(bulk(&c_sc) >= bulk(&c_r) - 1e6);
        assert!(bulk(&c_sc) <= bulk(&c_v) + 1e6);
        assert!(shear(&c_sc) >= shear(&c_r) - 1e6);
        assert!(shear(&c_sc) <= shear(&c_v) + 1e6);
        assert!(
            (bulk(&c_sc) - fcc_bulk).abs() < 0.02 * fcc_bulk,
            "SC bulk {} drifted from crystal bulk {fcc_bulk}",
            bulk(&c_sc)
        );
    }

    /// Bunge (Z–X–Z) Euler angles from a xorshift stream.
    fn random_bunge(seed: &mut u64) -> [[f64; 3]; 3] {
        let next = |seed: &mut u64| {
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;
            *seed as f64 / u64::MAX as f64
        };
        let phi1 = 2.0 * std::f64::consts::PI * next(seed);
        let phi = (1.0 - 2.0 * next(seed)).acos();
        let phi2 = 2.0 * std::f64::consts::PI * next(seed);
        let (s1, c1) = phi1.sin_cos();
        let (sp, cp) = phi.sin_cos();
        let (s2, c2) = phi2.sin_cos();
        [
            [c1 * c2 - s1 * cp * s2, s1 * c2 + c1 * cp * s2, sp * s2],
            [-c1 * s2 - s1 * cp * c2, -s1 * s2 + c1 * cp * c2, sp * c2],
            [s1 * sp, -c1 * sp, cp],
        ]
    }

    #[test]
    fn self_consistent_single_grain_recovers_input() {
        let c = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let rve = Rve::new(vec![RveGrain::new("g0", 1.0, identity(), c.clone())]);
        let c_sc = rve.homogenize(HomogenizationScheme::SelfConsistent);
        for i in 0..6 {
            for j in 0..6 {
                assert!((c_sc.data[i][j] - c.data[i][j]).abs() < 1e-3);
            }
        }
    }

    #[test]
    fn stats_counts_grains_and_volume_fraction() {
        let c = SymmetricFourthOrder::isotropic(100_000.0, 0.3);
        let rve = Rve::new(vec![
            RveGrain::new("a", 0.4, identity(), c.clone()),
            RveGrain::new("b", 0.6, identity(), c.clone()),
        ]);
        let s = rve.stats();
        assert_eq!(s.n_grains, 2);
        assert!((s.total_volume_fraction - 1.0).abs() < 1e-12);
        assert_eq!(s.unique_orientations, 1);
    }

    #[test]
    fn hill_mandel_voigt_uniform_strain_energy_consistency() {
        // Hill–Mandel macro-homogeneity for a uniform-strain
        // (Voigt) boundary: the macroscopic stress-energy
        // `<Sigma> : E` recovered from the homogenised stiffness
        // must equal the volume-averaged strain-energy
        // `<sigma : eps>` where each grain carries the same E.
        let c_a = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c_b = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let rve = Rve::new(vec![
            RveGrain::new("a", 0.3, identity(), c_a.clone()),
            RveGrain::new("b", 0.7, identity(), c_b.clone()),
        ]);
        let c_avg = rve.homogenize(HomogenizationScheme::Voigt);
        let e = Vec6::new(1.0e-3, -0.5e-3, 0.0, 0.0, 0.0, 0.5e-3);
        let sigma_avg = c_avg.contract(e);
        let sigma_a = c_a.contract(e);
        let sigma_b = c_b.contract(e);
        let macro_energy = sigma_avg.double_dot(e);
        let micro_energy = 0.3 * sigma_a.double_dot(e) + 0.7 * sigma_b.double_dot(e);
        let rel = (macro_energy - micro_energy).abs() / macro_energy.abs().max(1e-30);
        assert!(
            rel < 1.0e-12,
            "Hill–Mandel macro/micro energy mismatch: {rel:.3e}"
        );
    }

    #[test]
    fn hill_mandel_reuss_uniform_stress_energy_consistency() {
        // For Reuss (uniform-stress), the homogenised compliance
        // `S = <S^g>` yields macroscopic strain `E = S : Sigma`
        // for any applied stress Sigma.  The energy consistency
        // `<Sigma : E> = <Sigma : eps^g>` holds with all grains
        // seeing the same Sigma.
        let c_a = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c_b = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let phases = vec![(c_a.clone(), 0.3), (c_b.clone(), 0.7)];
        let s_avg = tpt_mat_homogenization::reuss(&phases);
        let s_a = c_a.compliance();
        let s_b = c_b.compliance();
        let sigma = Vec6::new(100.0, -50.0, 0.0, 0.0, 0.0, 25.0);
        let e_avg = apply_compliance(s_avg, sigma);
        let e_a = apply_compliance(s_a, sigma);
        let e_b = apply_compliance(s_b, sigma);
        let macro_energy = sigma.double_dot(e_avg);
        let micro_energy = sigma.double_dot(e_a.scale(0.3) + e_b.scale(0.7));
        let rel = (macro_energy - micro_energy).abs() / macro_energy.abs().max(1e-30);
        assert!(rel < 1.0e-9, "Reuss Hill–Mandel mismatch: {rel:.3e}");
    }
}

fn apply_compliance(s: [[f64; 6]; 6], v: Vec6) -> Vec6 {
    let mut out = [0.0_f64; 6];
    for i in 0..6 {
        let mut acc = 0.0;
        for j in 0..6 {
            acc += s[i][j] * v.data[j];
        }
        out[i] = acc;
    }
    Vec6::new(out[0], out[1], out[2], out[3], out[4], out[5])
}
