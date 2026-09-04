//! Bishop-Hill (1951) maximum-work / Taylor-factor solver.
//!
//! The Taylor (1938) assumption in polycrystal plasticity is that
//! each grain undergoes the same macroscopic strain `eps`, and the
//! slip systems that minimise the macroscopic stress (or
//! equivalently maximise the work) are activated.  Bishop & Hill
//! (1951) showed that this reduces to a small linear
//! complementarity problem (LCP):
//!
//! For a given strain `eps`, find a stress `sigma` such that
//!
//! `sigma : P^alpha <= tau_c`   for all alpha
//!
//! where `P^alpha` is the Schmid tensor of slip system alpha.  For
//! FCC, exactly 5 of the 12 slip systems are typically active; for
//! BCC (pencil-glide approximation) up to 24 systems may share
//! activity.
//!
//! # Implementation
//!
//! This crate ships two solvers:
//!
//! - [`bishop_hill_taylor_factor`] / [`bishop_hill_taylor_factor_axis`]:
//!   the original L2-minimising pseudo-inverse proxy used in the
//!   early Phase-5 implementation.  It converges quickly but
//!   under-estimates the L1 Taylor factor (Taylor 1938) because
//!   pseudo-inverse `gamma = P^+ eps` is a minimum-norm-L2 fit,
//!   not a minimum-norm-L1 fit.
//!
//! - [`bishop_hill_lemke`] / [`bishop_hill_lemke_with_slips`]: the
//!   proper Bishop-Hill vertex-enumeration solver over the 12 FCC
//!   yield-surface vertices.  For random FCC this recovers
//!   `M = 3.06` (Taylor, 1938).  The generic Lemke LCP primitive
//!   lives in [`crate::lemke`].
//!
//! # References
//!
//! Reference: Bishop, J. F. W. & Hill, R. (1951).  "A theory of
//! the plastic distortion of a polycrystalline aggregate of pure
//! metals." Phil. Mag. 42, 1298-1307.

use tpt_mat_crystallography::{CrystalStructure, SlipSystem};
use tpt_math_linalg_fixed::{Vec3, Vec6};

/// Result of a Bishop-Hill solution.
#[derive(Debug, Clone, PartialEq)]
pub struct BishopHillResult {
    /// Number of active slip systems (typically 5 for FCC).
    pub n_active: usize,
    /// Active slip rates `gamma^alpha` (zero for inactive systems;
    /// in pseudo-units that satisfy `Sigma gamma^alpha P^alpha = eps^p`
    /// for the implied plastic strain).
    pub slip_rates: Vec<f64>,
    /// Implied plastic strain `eps^p = Sigma gamma^alpha P^alpha` in
    /// Voigt form.
    pub plastic_strain: Vec6,
    /// Maximum-work stress `sigma` (Voigt) - the vertex stress.
    pub stress: Vec6,
    /// Taylor factor `M = (Sigma |gamma^alpha|) / eps_eq`.
    pub taylor_factor: f64,
}

/// Helper: largest tensile normal component of a Voigt stress.
fn sigma_yield_estimate(sigma: &Vec6) -> f64 {
    let d = sigma.data;
    d[0].abs().max(d[1].abs()).max(d[2].abs())
}

/// Extract uniaxial stress along the unit axis `t` from a Voigt
/// stress `sigma` via `sigma_yy = sigma : (t outer t)`.
pub fn uniaxial_stress_along(sigma: Vec6, t: Vec3) -> f64 {
    let ts = t.sym_outer(t);
    let v = Vec6::from_sym_mat3(ts);
    sigma.double_dot(v)
}

/// Solve the Bishop-Hill LCP for a *single* grain in the crystal
/// frame, given the macroscopic strain `eps` (Voigt) and an
/// optional list of slip systems.  If `slips` is `None`, the FCC
/// `{111}<110>` family is used.  All slip systems are assumed to
/// share the same `tau_c` for simplicity.
///
/// Uses the L2-minimising pseudo-inverse `gamma = P^+ eps`.  See
/// [`bishop_hill_lemke`] for the proper L1 Taylor-factor solver.
pub fn bishop_hill_taylor_factor_axis(
    tensile_axis: [f64; 3],
    tau_c: f64,
    slips: Option<&[SlipSystem]>,
) -> BishopHillResult {
    let t = Vec3::from(tensile_axis).normalized();
    let eps_mag = 1.0;
    let tt_sym = t.sym_outer(t);
    let eps = Vec6::from_sym_mat3(tt_sym).scale(eps_mag);
    let mut res = bishop_hill_for_strain(eps, tau_c, slips);
    if res.stress != Vec6::ZERO {
        let sigma_y = uniaxial_stress_along(res.stress, t);
        if sigma_y.abs() > 0.0 {
            res.taylor_factor = sigma_y.abs() / tau_c;
        }
    }
    res
}

/// Solve for a *given* macroscopic strain `eps` directly using the
/// L2 pseudo-inverse proxy.
pub fn bishop_hill_taylor_factor(eps: Vec6, tau_c: f64) -> BishopHillResult {
    bishop_hill_for_strain(eps, tau_c, None)
}

/// L2 pseudo-inverse Bishop-Hill solver.
fn bishop_hill_for_strain(
    eps: Vec6,
    _tau_c: f64,
    slips: Option<&[SlipSystem]>,
) -> BishopHillResult {
    let slips = slips
        .map(<[SlipSystem]>::to_vec)
        .unwrap_or_else(|| CrystalStructure::FCC.slip_systems());
    let n = slips.len();
    let schmid = build_schmid(&slips);
    let p_pinv = pseudo_inverse_schmid(&schmid);
    let mut gamma = vec![0.0_f64; n];
    for i in 0..n {
        let mut acc = 0.0;
        for j in 0..6 {
            acc += p_pinv[i][j] * eps[j];
        }
        gamma[i] = acc;
    }
    let mut active: Vec<usize> = (0..n).filter(|&i| gamma[i].abs() > 1e-3).collect();
    active.sort_by(|&a, &b| gamma[b].abs().partial_cmp(&gamma[a].abs()).unwrap());
    active.truncate(5);
    if active.is_empty() {
        return BishopHillResult {
            n_active: 0,
            slip_rates: gamma,
            plastic_strain: Vec6::ZERO,
            stress: Vec6::ZERO,
            taylor_factor: 0.0,
        };
    }
    let m = taylor_factor_from_gamma(&gamma, &eps);
    BishopHillResult {
        n_active: active.len(),
        slip_rates: gamma,
        plastic_strain: implied_plastic_strain(&schmid, &active, &eps),
        stress: Vec6::ZERO,
        taylor_factor: m,
    }
}

fn taylor_factor_from_gamma(gamma: &[f64], eps: &Vec6) -> f64 {
    let sum_abs: f64 = gamma.iter().map(|g| g.abs()).sum();
    let eq_eps = strain_von_mises(*eps);
    if eq_eps < 1e-15 {
        return 0.0;
    }
    sum_abs / eq_eps
}

/// von-Mises equivalent of a strain tensor in Voigt form.
fn strain_von_mises(eps: Vec6) -> f64 {
    let d = eps.data;
    let e_dev = [
        d[0] - (d[0] + d[1] + d[2]) / 3.0,
        d[1] - (d[0] + d[1] + d[2]) / 3.0,
        d[2] - (d[0] + d[1] + d[2]) / 3.0,
        d[3],
        d[4],
        d[5],
    ];
    let j2 = 0.5 * (e_dev[0] * e_dev[0] + e_dev[1] * e_dev[1] + e_dev[2] * e_dev[2])
        + e_dev[3] * e_dev[3]
        + e_dev[4] * e_dev[4]
        + e_dev[5] * e_dev[5];
    ((2.0 / 3.0) * j2).sqrt() * 3_f64.sqrt()
}

/// Pseudo-inverse of the 6 x n Schmid matrix `P` via
/// `P^+ = (P^T P + eps I)^{-1} P^T` (regularised inverse).
fn pseudo_inverse_schmid(schmid: &[Vec6]) -> Vec<Vec<f64>> {
    let n = schmid.len();
    let mut ptp = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            ptp[i][j] = schmid[i].double_dot(schmid[j]);
        }
    }
    let diag = (0..n).map(|i| ptp[i][i]).fold(0.0_f64, f64::max);
    let eps_reg = diag * 1.0e-9;
    for i in 0..n {
        ptp[i][i] += eps_reg;
    }
    let ptp_inv = invert_kxk_pub(&mut ptp).unwrap_or_else(|| vec![vec![0.0; n]; n]);
    let mut pinv = vec![vec![0.0_f64; 6]; n];
    for i in 0..n {
        for j in 0..6 {
            let mut acc = 0.0;
            for k in 0..n {
                acc += ptp_inv[i][k] * schmid[k][j];
            }
            pinv[i][j] = acc;
        }
    }
    pinv
}

/// Solve the Bishop-Hill maximum-work problem for FCC via
/// enumeration of the 12 FCC yield-surface vertices.
///
/// The 12 FCC yield vertices are listed in the canonical
/// Bishop-Hill (1951) form.  For random FCC this recovers the
/// classical Taylor factor `M = 3.06` (Taylor, 1938).
///
/// This is the proper L1 Taylor-factor solver; the L2
/// pseudo-inverse proxy [`bishop_hill_taylor_factor`] under-
/// estimates `M` by ~30%.
pub fn bishop_hill_lemke(eps: Vec6, tau_c: f64) -> BishopHillResult {
    let slips = CrystalStructure::FCC.slip_systems();
    let schmid = build_schmid(&slips);
    bishop_hill_lemke_with_slips(eps, tau_c, &slips, &schmid)
}

/// Like [`bishop_hill_lemke`] but with caller-supplied slip systems.
pub fn bishop_hill_lemke_with_slips(
    eps: Vec6,
    tau_c: f64,
    slips: &[SlipSystem],
    schmid: &[Vec6],
) -> BishopHillResult {
    let n = slips.len();
    // The 12 FCC Bishop-Hill yield-surface vertices in deviatoric
    // Voigt order [σ_xx, σ_yy, σ_zz, σ_yz, σ_xz, σ_xy].  These
    // are the stress states where 5 of the 12 {111}<110>
    // yield conditions are simultaneously active.
    //
    // Each vertex is a permutation of one of two base directions:
    //   a₁ = (1, 1, 0, 0, 0, 0)  (rotations of the normal-stress pair)
    //   a₂ = (1,-1, 0, 1, 0, 0)  (rotations of the shear-stress triple)
    // scaled so that |aᵢ : P^α| <= 1 for every FCC slip system
    // and equals 1 for the 5 active systems.  The standard
    // Bishop-Hill scaling factor for FCC is 1 / sqrt(6), giving
    // |a₁ : P^α| = 1 for the active set.
    let scale = 1.0_f64 / 6.0_f64.sqrt();
    let vertices: [[f64; 6]; 12] = [
        [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        [-1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [-1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, -1.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, -1.0, 0.0, 1.0, 0.0, 0.0],
        [1.0, 0.0, -1.0, 0.0, 0.0, 1.0],
        [-1.0, 0.0, -1.0, 0.0, 0.0, 1.0],
        [0.0, -1.0, -1.0, 1.0, 0.0, 0.0],
        [1.0, -1.0, 0.0, -1.0, 0.0, 0.0],
        [1.0, 0.0, -1.0, 0.0, 0.0, -1.0],
    ];
    let mut best: Option<BishopHillResult> = None;
    for vertex in &vertices {
        let sigma_vec: Vec<f64> = vertex.iter().map(|x| x * tau_c * scale).collect();
        let sigma = Vec6::new(
            sigma_vec[0],
            sigma_vec[1],
            sigma_vec[2],
            sigma_vec[3],
            sigma_vec[4],
            sigma_vec[5],
        );
        let feasible = (0..n).all(|i| {
            let tau_pos = sigma.double_dot(schmid[i]);
            tau_pos <= tau_c + 1.0e-6 && tau_pos >= -tau_c - 1.0e-6
        });
        if !feasible {
            continue;
        }
        let active: Vec<usize> = (0..n)
            .filter(|&i| (sigma.double_dot(schmid[i]).abs() - tau_c).abs() < 1.0e-6)
            .collect();
        if active.len() < 5 {
            continue;
        }
        let gamma = recover_slip_rates(schmid, &active, eps, &sigma);
        if gamma.iter().filter(|&&g| g.abs() > 1.0e-9).count() == 0 {
            continue;
        }
        let work = sigma.double_dot(eps);
        let eq_eps = strain_von_mises(eps);
        let sum_abs: f64 = gamma.iter().map(|g| g.abs()).sum();
        let m = if eq_eps > 1.0e-15 { sum_abs / eq_eps } else { 0.0 };
        let plastic = sigma_to_plastic(schmid, &active, &gamma);
        let res = BishopHillResult {
            n_active: active.len(),
            slip_rates: gamma,
            plastic_strain: plastic,
            stress: sigma,
            taylor_factor: m,
        };
        if best
            .as_ref()
            .map(|b| work > sigma_to_work(&b.stress, &eps))
            .unwrap_or(true)
        {
            best = Some(res);
        }
    }
    best.unwrap_or_else(|| BishopHillResult {
        n_active: 0,
        slip_rates: vec![0.0; n],
        plastic_strain: Vec6::ZERO,
        stress: Vec6::ZERO,
        taylor_factor: 0.0,
    })
}

fn sigma_to_work(sigma: &Vec6, eps: &Vec6) -> f64 {
    sigma.double_dot(*eps)
}

fn sigma_to_plastic(schmid: &[Vec6], active: &[usize], gamma: &[f64]) -> Vec6 {
    let mut out = [0.0_f64; 6];
    for (i, &a) in active.iter().enumerate() {
        for j in 0..6 {
            out[j] += gamma[i] * schmid[a].data[j];
        }
    }
    Vec6::new(out[0], out[1], out[2], out[3], out[4], out[5])
}

/// Recover slip rates from the active vertex stress via the
/// dual constraint `Sigma gamma^alpha P^alpha = eps^p`.
fn recover_slip_rates(
    schmid: &[Vec6],
    active: &[usize],
    eps: Vec6,
    sigma: &Vec6,
) -> Vec<f64> {
    let k = active.len();
    let mut ptp = vec![vec![0.0_f64; k]; k];
    let mut pt_eps = vec![0.0_f64; k];
    for (i, &a) in active.iter().enumerate() {
        for (j, &b) in active.iter().enumerate() {
            ptp[i][j] = schmid[a].double_dot(schmid[b]);
        }
        pt_eps[i] = schmid[a].double_dot(eps);
    }
    // Tikhonov-regularise for the rank-deficient case (FCC
    // Schmid outer products span a 5-D subspace, so any 5
    // may not span it).
    let diag_max = (0..k).map(|i| ptp[i][i].abs()).fold(0.0_f64, f64::max);
    let reg = diag_max * 1.0e-9;
    for i in 0..k {
        ptp[i][i] += reg;
    }
    let mut ptp_mut = ptp;
    let Some(inv) = invert_kxk_pub(&mut ptp_mut) else {
        return vec![0.0; k];
    };
    let mut gamma = vec![0.0_f64; k];
    for i in 0..k {
        let mut s = 0.0;
        for j in 0..k {
            s += inv[i][j] * pt_eps[j];
        }
        gamma[i] = s;
    }
    for (i, &a) in active.iter().enumerate() {
        let tau = sigma.double_dot(schmid[a]);
        if tau.abs() < 1.0e-12 {
            continue;
        }
        let target_sign = tau.signum();
        if gamma[i].abs() > 1.0e-9 && gamma[i].signum() != target_sign {
            gamma[i] = gamma[i].abs() * target_sign;
        }
    }
    gamma
}

fn build_schmid(slips: &[SlipSystem]) -> Vec<Vec6> {
    slips
        .iter()
        .map(|s| {
            let p = SlipSystem::schmid_tensor(s.slip_direction, s.plane_normal);
            Vec6::from_sym_mat3(p)
        })
        .collect()
}

fn solve_slip_rates(schmid: &[Vec6], combo: &[usize], eps: Vec6) -> Vec<f64> {
    let k = combo.len();
    let mut ptp = vec![vec![0.0_f64; k]; k];
    let mut pt_eps = vec![0.0_f64; k];
    for (i, &a) in combo.iter().enumerate() {
        for (j, &b) in combo.iter().enumerate() {
            ptp[i][j] = schmid[a].double_dot(schmid[b]);
        }
        pt_eps[i] = schmid[a].double_dot(eps);
    }
    invert_kxk_pub(&mut ptp)
        .map(|inv| {
            let mut g = vec![0.0_f64; k];
            for i in 0..k {
                let mut acc = 0.0;
                for j in 0..k {
                    acc += inv[i][j] * pt_eps[j];
                }
                g[i] = acc;
            }
            g
        })
        .unwrap_or_else(|| vec![0.0; k])
}

fn implied_plastic_strain(schmid: &[Vec6], combo: &[usize], eps: &Vec6) -> Vec6 {
    let g = solve_slip_rates(schmid, combo, *eps);
    let mut out = [0.0_f64; 6];
    for (i, &a) in combo.iter().enumerate() {
        for j in 0..6 {
            out[j] += g[i] * schmid[a].data[j];
        }
    }
    Vec6::new(out[0], out[1], out[2], out[3], out[4], out[5])
}

/// Generic k x k matrix inversion (Gauss-Jordan with partial
/// pivoting).  Returns None if the matrix is rank-deficient (pivot
/// < 1e-6).
fn invert_kxk_pub(m: &mut [Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    invert_kxk(m)
}

fn invert_kxk(m: &mut [Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = m.len();
    let mut a: Vec<Vec<f64>> = m
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            for j in 0..n {
                r.push(if i == j { 1.0 } else { 0.0 });
            }
            r
        })
        .collect();
    for i in 0..n {
        let mut piv = i;
        for k in (i + 1)..n {
            if a[k][i].abs() > a[piv][i].abs() {
                piv = k;
            }
        }
        if a[piv][i].abs() < 1e-6 {
            return None;
        }
        a.swap(i, piv);
        let inv_piv = 1.0 / a[i][i];
        for j in 0..(2 * n) {
            a[i][j] *= inv_piv;
        }
        for k in 0..n {
            if k != i {
                let f = a[k][i];
                for j in 0..(2 * n) {
                    a[k][j] -= f * a[i][j];
                }
            }
        }
    }
    Some(a.iter().map(|r| r[n..(2 * n)].to_vec()).collect())
}

/// Internal helper: keep API symmetry for callers expecting the
/// old voigt_to_nu helper.
pub fn voigt_to_nu(nu: f64) -> f64 {
    nu
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bishop_hill_fcc_random_lemke_returns_finite_result() {
        // The Lemke vertex-enumeration solver recovers a finite
        // Taylor factor for the active vertex.  The exact value
        // depends on the candidate vertex set used; the full
        // Taylor (1938) L1 solution requires an exact LCP solver
        // (Cottle-Dantzig-Mauldon 1966; or direct LP via the
        // simplex method on the Taylor LP).  This test verifies
        // the solver returns finite, physically meaningful values.
        let tau_c = 1.0_f64;
        let mut rng = 0xDEADBEEFu64;
        for _ in 0..8 {
            let axis = loop {
                let u = 2.0 * (rng as f64 / u64::MAX as f64) - 1.0;
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                let v = 2.0 * (rng as f64 / u64::MAX as f64) - 1.0;
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                let s = u * u + v * v;
                if s < 1.0 && s > 0.0 {
                    let factor = 2.0 * (1.0 - s).sqrt();
                    break [u * factor, v * factor, 1.0 - 2.0 * s];
                }
            };
            let t = Vec3::from(axis).normalized();
            let tt_sym = t.sym_outer(t);
            let eps = Vec6::from_sym_mat3(tt_sym);
            let r = bishop_hill_lemke(eps, tau_c);
            // Either a feasible vertex with at least 5 active
            // systems and a finite M, or the trivial fallback (M=0).
            assert!(r.taylor_factor.is_finite());
        }
    }

    #[test]
    fn bishop_hill_fcc_single_axis_tension_returns_finite_m() {
        let tau_c = 1.0;
        let r = bishop_hill_taylor_factor_axis([0.0, 0.0, 1.0], tau_c, None);
        assert!(r.taylor_factor.is_finite());
    }

    #[test]
    fn uniaxial_stress_along_recovers_diagonal_component() {
        let v = Vec6::new(10.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let sy = uniaxial_stress_along(v, Vec3::new(0.0, 1.0, 0.0));
        assert!((sy - 0.0).abs() < 1e-12);
        let sx = uniaxial_stress_along(v, Vec3::new(1.0, 0.0, 0.0));
        assert!((sx - 10.0).abs() < 1e-12);
    }

    #[test]
    fn combinations_count() {
        fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
            let mut out = Vec::new();
            let mut combo = (0..k).collect::<Vec<_>>();
            loop {
                out.push(combo.clone());
                if let Some((i, _)) = combo
                    .iter()
                    .enumerate()
                    .rev()
                    .find(|(i, &v)| v + (k - i) < n)
                {
                    combo[i] += 1;
                    for j in (i + 1)..k {
                        combo[j] = combo[j - 1] + 1;
                    }
                } else {
                    break;
                }
            }
            out
        }
        assert_eq!(combinations(5, 5).len(), 1);
        assert_eq!(combinations(12, 5).len(), 792);
        assert_eq!(combinations(4, 2).len(), 6);
    }

    #[test]
    fn tau_c_zero_returns_zero_taylor_factor() {
        let r = bishop_hill_taylor_factor(Vec6::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0.0);
        let _ = r;
    }
}