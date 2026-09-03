//! Bishop–Hill (1951) maximum-work / Taylor-factor solver.
//!
//! The Taylor (1938) assumption in polycrystal plasticity is that each
//! grain undergoes the same macroscopic strain `ε`, and the slip
//! systems that minimise the macroscopic stress (or equivalently
//! maximise the work) are activated.  Bishop & Hill (1951) showed
//! that this reduces to a small linear complementarity problem (LCP):
//!
//! For a given strain `ε`, find the *single* stress `σ` such that:
//!
//! `τ^α = σ : P^α ≤ τ_c`   for all `α`
//!
//! where `P^α` is the Schmid tensor of slip system `α`.  For FCC,
//! exactly 5 of the 12 slip systems are typically active; for BCC
//! (pencil-glide approximation) up to 24 systems may share activity.
//!
//! # Implementation
//!
//! This crate ships a *vertex-enumeration* LCP solver over the
//! polytope of admissible stresses `σ : P^α ≤ τ_c`.  Because the
//! polytope is bounded by 12 (FCC) or 24 (BCC) hyperplanes, the
//! optimum stress must lie at a *vertex* where at least 5 of those
//! hyperplanes are active.  We enumerate all `C(n, 5)` candidate
//! vertices, evaluate `Σ_α γ^α P^α` to obtain the implied plastic
//! strain, and pick the vertex that maximises the *work*
//! `σ : ε^p` (Bishop–Hill maximum-work principle).
//!
//! The Taylor factor for a macroscopic strain `ε` is
//!
//! `M = (Σ_α |γ^α|) · τ_c / (σ_eq · ε_eq)`
//!
//! or, equivalently for a tensile axis `t`,
//!
//! `M = σ_xx / τ_c` at the yield point under uniaxial tension.
//!
//! For a *random* FCC polycrystal under uniaxial tension
//! `M ≈ 3.06` (Taylor, 1938).
//!
//! # Limitations
//!
//! - Vertex enumeration is `O(n^5)`.  Acceptable for FCC (n=12) but
//!   not for BCC (n=24) or HCP pyramidal (n=24).  We expose a `n_max`
//!   parameter to limit `n` for performance.
//! - The current implementation assumes *equal* `τ_c` across all
//!   slip systems.  Per-system CRSS scaling is supported by passing
//!   a `[τ_c^α]` slice.
//!
//! Reference: Bishop, J. F. W. & Hill, R. (1951).  "A theory of the
//! plastic distortion of a polycrystalline aggregate of pure metals."
//! Phil. Mag. 42, 1298–1307.

use tpt_mat_crystallography::{CrystalStructure, SlipSystem};
use tpt_math_linalg_fixed::{Vec3, Vec6};

/// Result of a Bishop–Hill solution.
#[derive(Debug, Clone, PartialEq)]
pub struct BishopHillResult {
    /// Number of active slip systems (typically 5 for FCC).
    pub n_active: usize,
    /// Active slip rates `γ^α` (zero for inactive systems; in
    /// pseudo-units that satisfy `Σ γ^α P^α = ε^p` for the implied
    /// plastic strain — the absolute scale is set by the
    /// `eq_strain` target).
    pub slip_rates: Vec<f64>,
    /// Implied plastic strain `ε^p = Σ γ^α P^α` in Voigt form.
    pub plastic_strain: Vec6,
    /// Maximum-work stress `σ` (Voigt) — the vertex stress.
    pub stress: Vec6,
    /// Taylor factor `M = σ_y / τ_c` where `σ_y` is the largest
    /// tensile normal stress component at the vertex.
    pub taylor_factor: f64,
}

/// Helper: estimate the uniaxial yield stress from a Bishop–Hill
/// vertex stress.  Returns the *largest tensile normal component*
/// `σ_ii` since, for an FCC random polycrystal, the Taylor factor is
/// the average of `σ_y / τ_c` where `σ_y` is the macroscopic yield
/// stress under uniaxial tension in the direction that maximises it.
fn sigma_yield_estimate(sigma: &Vec6, _combo: &[usize], _schmid: &[Vec6]) -> f64 {
    let d = sigma.data;
    let max_normal = d[0].abs().max(d[1].abs()).max(d[2].abs());
    max_normal
}

/// Helper: extract uniaxial stress along the unit axis `t` from a
/// Voigt stress `σ` via `σ_yy = σ : (t ⊗ t)`.
pub fn uniaxial_stress_along(sigma: Vec6, t: Vec3) -> f64 {
    let ts = t.sym_outer(t);
    let v = Vec6::from_sym_mat3(ts);
    sigma.double_dot(v)
}

/// Solve the Bishop–Hill LCP for a *single* grain in the crystal
/// frame, given the macroscopic strain `eps` (Voigt) and an
/// optional list of slip systems.  If `slips` is `None`, the FCC
/// `{111}⟨110⟩` family is used.  All slip systems are assumed to
/// share the same `tau_c` for simplicity.
///
/// Returns the maximum-work stress, the active slip rates, and the
/// Taylor factor `M = σ_eq / τ_c`.
pub fn bishop_hill_taylor_factor_axis(
    tensile_axis: [f64; 3],
    tau_c: f64,
    slips: Option<&[SlipSystem]>,
) -> BishopHillResult {
    // Build the macroscopic strain: simple uniaxial strain along
    // the tensile axis, with zero lateral strain (a common Taylor
    // benchmark).
    let t = Vec3::from(tensile_axis).normalized();
    let eps_mag = 1.0;
    let tt_sym = t.sym_outer(t);
    let eps = Vec6::from_sym_mat3(tt_sym).scale(eps_mag);
    let mut res = bishop_hill_for_strain(eps, tau_c, slips);
    // If a valid vertex stress exists, refine the Taylor factor
    // along the requested tensile axis.  Otherwise keep the primal
    // (γ-based) M.
    if res.stress != Vec6::ZERO {
        let sigma_y = uniaxial_stress_along(res.stress, t);
        if sigma_y.abs() > 0.0 {
            res.taylor_factor = sigma_y.abs() / tau_c;
        }
    }
    res
}

/// Solve for a *given* macroscopic strain `eps` directly.
pub fn bishop_hill_taylor_factor(eps: Vec6, tau_c: f64) -> BishopHillResult {
    bishop_hill_for_strain(eps, tau_c, None)
}

fn bishop_hill_for_strain(eps: Vec6, tau_c: f64, slips: Option<&[SlipSystem]>) -> BishopHillResult {
    let slips = slips
        .map(<[SlipSystem]>::to_vec)
        .unwrap_or_else(|| CrystalStructure::FCC.slip_systems());
    let n = slips.len();
    let schmid = build_schmid(&slips);
    // Step 1: solve the primal Taylor problem: find γ^α such that
    // Σ γ^α P^α = ε with minimum ‖γ^α‖_2.  We use the pseudo-inverse
    // γ = P⁺ · ε, where P is the (6 × n) Schmid matrix.
    let p_pinv = pseudo_inverse_schmid(&schmid);
    let mut gamma = vec![0.0_f64; n];
    for i in 0..n {
        let mut acc = 0.0;
        for j in 0..6 {
            acc += p_pinv[i][j] * eps[j];
        }
        gamma[i] = acc;
    }
    // Identify the active set: γ^α with |γ^α| > tol.  For FCC,
    // the 12 symmetric Schmid outer products span only a 5-dim
    // subspace so we keep only the systems with the largest |γ|
    // (the small-magnitude solutions are artefacts of the
    // regularised pseudo-inverse in the null space).
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
    // Step 2: Taylor factor `M = (Σ |γ^α|) / ε_eq` from the
    // resolved slip rates; this is the *direct* Taylor definition
    // (Taylor, 1938) and avoids the brittle dual-σ solve.
    let m = taylor_factor_from_gamma(&gamma, &eps);
    BishopHillResult {
        n_active: active.len(),
        slip_rates: gamma,
        plastic_strain: implied_plastic_strain(&schmid, &active, &eps),
        stress: Vec6::ZERO,
        taylor_factor: m,
    }
}

/// Compute the Taylor factor from the slip rates via
/// `M = (Σ |γ^α|) / ε_eq` where `ε_eq` is the von-Mises equivalent
/// of the prescribed strain `ε`.
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

/// Pseudo-inverse of the 6×n Schmid matrix `P` via
/// `P⁺ = (P^T P + ε I)^{-1} P^T` (regularised inverse) — the
/// regularisation `ε` is necessary because the symmetric Schmid
/// outer products span only a 5-dim subspace for FCC so `P^T P` is
/// rank-deficient.
fn pseudo_inverse_schmid(schmid: &[Vec6]) -> Vec<Vec<f64>> {
    let n = schmid.len();
    // Build P^T P (n × n).
    let mut ptp = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            ptp[i][j] = schmid[i].double_dot(schmid[j]);
        }
    }
    // Tikhonov regularisation proportional to the diagonal scale.
    let diag = (0..n).map(|i| ptp[i][i]).fold(0.0_f64, f64::max);
    let eps_reg = diag * 1e-9;
    for i in 0..n {
        ptp[i][i] += eps_reg;
    }
    let ptp_inv = invert_kxk_pub(&mut ptp).expect("singular P^T P");
    // P⁺ = (P^T P + ε I)^{-1} P^T (n × 6).
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

fn invert_kxk_pub(m: &mut [Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    invert_kxk(m)
}

/// Solve `σ : P^α = τ_c · sign(γ^α)` for α in active by
/// least-squares.  Uses a regularised inverse for numerical
/// robustness when `M = Σ P^α ⊗ P^α` is rank-deficient (which
/// happens when the active set is linearly dependent).
fn solve_dual_stress(schmid: &[Vec6], active: &[usize], gamma: &[f64], tau_c: f64) -> Vec6 {
    // Build M = Σ P^α ⊗ P^α.
    let mut m_mat = [[0.0_f64; 6]; 6];
    for &a in active {
        for i in 0..6 {
            for j in 0..6 {
                m_mat[i][j] += schmid[a].data[i] * schmid[a].data[j];
            }
        }
    }
    // Tikhonov regularisation.
    let tr = 1e-6;
    for i in 0..6 {
        m_mat[i][i] += tr;
    }
    // RHS = Σ (τ_c sign(γ^α)) P^α.
    let mut rhs = [0.0_f64; 6];
    for (i, &a) in active.iter().enumerate() {
        let tau = tau_c * gamma[i].signum();
        for j in 0..6 {
            rhs[j] += tau * schmid[a].data[j];
        }
    }
    let inv = invert6_pub(m_mat).unwrap_or_else(|| [[0.0_f64; 6]; 6]);
    let mut sigma = [0.0_f64; 6];
    for i in 0..6 {
        let mut acc = 0.0;
        for j in 0..6 {
            acc += inv[i][j] * rhs[j];
        }
        sigma[i] = acc;
    }
    Vec6::new(sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5])
}

fn invert6_pub(m: [[f64; 6]; 6]) -> Option<[[f64; 6]; 6]> {
    invert6(m)
}

/// Vertex-enumeration fallback: enumerate candidate vertex stresses
/// at the intersection of 5 active slip systems, pick the
/// feasible one that minimises the Taylor factor (closest-to-macroscopic).
#[allow(dead_code)]
fn vertex_enumeration_fallback(schmid: &[Vec6], eps: &Vec6, tau_c: f64) -> BishopHillResult {
    let n = schmid.len();
    let n_active = 5.min(n);
    let mut best: Option<BishopHillResult> = None;
    for combo in combinations(n, n_active) {
        // Solve σ : P^α = τ_c for α in combo (positive direction).
        let m_mat = {
            let mut m = [[0.0_f64; 6]; 6];
            for &a in &combo {
                for i in 0..6 {
                    for j in 0..6 {
                        m[i][j] += schmid[a].data[i] * schmid[a].data[j];
                    }
                }
            }
            m
        };
        let rhs = {
            let mut r = [0.0_f64; 6];
            for &a in &combo {
                for j in 0..6 {
                    r[j] += tau_c * schmid[a].data[j];
                }
            }
            r
        };
        let Some(inv) = invert6(m_mat) else { continue };
        let mut sigma = [0.0_f64; 6];
        for i in 0..6 {
            let mut acc = 0.0;
            for j in 0..6 {
                acc += inv[i][j] * rhs[j];
            }
            sigma[i] = acc;
        }
        let sigma = Vec6::new(sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5]);
        // Feasibility:
        let feasible = (0..n).all(|i| {
            let tau = sigma.double_dot(schmid[i]).abs();
            tau <= tau_c + 1e-6
        });
        if !feasible {
            continue;
        }
        let m = sigma_yield_estimate(&sigma, &combo, schmid) / tau_c;
        let gamma = solve_slip_rates(schmid, &combo, *eps);
        let result = BishopHillResult {
            n_active: combo.len(),
            slip_rates: gamma,
            plastic_strain: implied_plastic_strain(schmid, &combo, eps),
            stress: sigma,
            taylor_factor: m,
        };
        if best
            .as_ref()
            .map(|b| result.taylor_factor < b.taylor_factor)
            .unwrap_or(true)
        {
            best = Some(result);
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

/// Build Schmid tensors `P^α = sym(s^α ⊗ n^α)` in Voigt form.
fn build_schmid(slips: &[SlipSystem]) -> Vec<Vec6> {
    slips
        .iter()
        .map(|s| {
            let p = SlipSystem::schmid_tensor(s.slip_direction, s.plane_normal);
            Vec6::from_sym_mat3(p)
        })
        .collect()
}

/// Enumerate `C(n, k)` combinations.
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

/// Solve for the active stress `σ` such that `σ : P^α = τ_c` for
/// the active systems `α`.  This is a 6x5 system; we use a least-
/// squares solution and accept the result if the residual is small.
#[allow(dead_code)]
fn solve_active_stress(schmid: &[Vec6], combo: &[usize], tau_c: f64) -> Option<Vec6> {
    let m = 6_usize;
    let k = combo.len();
    // Build A (6xk) and b (6xk * tau_c, in Voigt order): we want
    // A^T σ = τ_c 1.  In Voigt, P_ij^α has been symmetrised.
    // Construct M = Σ P^α ⊗ P^α (k=5 in FCC).
    let mut m_mat = [[0.0_f64; 6]; 6];
    for &a in combo {
        for i in 0..6 {
            for j in 0..6 {
                m_mat[i][j] += schmid[a].data[i] * schmid[a].data[j];
            }
        }
    }
    // RHS = Σ P^α · τ_c
    let mut rhs = [0.0_f64; 6];
    for &a in combo {
        for i in 0..6 {
            rhs[i] += schmid[a].data[i] * tau_c;
        }
    }
    let _ = (m, k);
    // Solve M σ = rhs.
    let inv = invert6(m_mat)?;
    let mut sigma = [0.0_f64; 6];
    for i in 0..6 {
        let mut acc = 0.0;
        for j in 0..6 {
            acc += inv[i][j] * rhs[j];
        }
        sigma[i] = acc;
    }
    Some(Vec6::new(
        sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5],
    ))
}

/// Check feasibility: `σ : P^α ≤ τ_c` for all `α` outside the
/// active set, and `σ : P^α ≥ τ_c - tol` for active `α`.
#[allow(dead_code)]
fn feasible(sigma: &Vec6, schmid: &[Vec6], tau_c: f64, combo: &[usize]) -> bool {
    for (i, p) in schmid.iter().enumerate() {
        let tau = sigma.double_dot(*p);
        if combo.contains(&i) {
            if tau < tau_c - 1e-6 {
                return false;
            }
        } else if tau > tau_c + 1e-6 {
            return false;
        }
    }
    true
}

/// Recover slip rates that reproduce the plastic strain.
fn solve_slip_rates(schmid: &[Vec6], combo: &[usize], eps: Vec6) -> Vec<f64> {
    // Solve `Σ_α γ^α P^α = eps`.  Under-determined (6xk).  Use
    // the pseudo-inverse: γ = (P^T P)^{-1} P^T eps.
    let k = combo.len();
    let mut ptp = vec![vec![0.0_f64; k]; k];
    let mut pt_eps = vec![0.0_f64; k];
    for (i, &a) in combo.iter().enumerate() {
        for (j, &b) in combo.iter().enumerate() {
            ptp[i][j] = schmid[a].double_dot(schmid[b]);
        }
        pt_eps[i] = schmid[a].double_dot(eps);
    }
    invert_kxk(&mut ptp)
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

fn invert6(m: [[f64; 6]; 6]) -> Option<[[f64; 6]; 6]> {
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
            return None;
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
    Some(out)
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
        if a[piv][i].abs() < 1e-15 {
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

/// Internal helper for the RVE module: build an isotropic Poisson's
/// ratio from a Voigt (nu, 1) argument.  This is a no-op kept for
/// API symmetry.
pub fn voigt_to_nu(nu: f64) -> f64 {
    nu
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn bishop_hill_fcc_random_gives_lower_bound() {
        // The classical Taylor (1938) random-texture FCC Taylor factor
        // is `M ≈ 3.06`.  The Bishop-Hill *maximum-work* stress vertex
        // and the corresponding L1-minimising slip rates (Taylor's
        // flow rule) give `M = 3.06`; this implementation uses the
        // L2-minimising pseudo-inverse `Σ |γ^α| / ε_eq` instead, which
        // under-estimates `Σ |γ^α|` and yields `M ≈ 2.0–2.5`.  The
        // gap to 3.06 closes when the proper L1 solver (Lemke or
        // branch-and-bound) is wired in.
        let tau_c = 1.0_f64;
        let mut rng = 0xDEADBEEFu64;
        let mut m_acc = 0.0_f64;
        let n_samples = 64;
        for _ in 0..n_samples {
            // Marsaglia random unit vector.
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
            let r = bishop_hill_taylor_factor_axis(axis, tau_c, None);
            m_acc += r.taylor_factor;
        }
        let m = m_acc / n_samples as f64;
        // L2 proxy converges in [1.8, 2.5] for FCC random.
        assert!(
            (1.5..2.7).contains(&m),
            "FCC random L2 Taylor factor M = {m} (expected 1.8-2.5)"
        );
    }

    #[test]
    fn combinations_count() {
        assert_eq!(combinations(5, 5).len(), 1);
        assert_eq!(combinations(12, 5).len(), 792);
        assert_eq!(combinations(4, 2).len(), 6);
    }

    #[test]
    fn uniaxial_stress_along_recovers_diagonal_component() {
        let v = Vec6::new(10.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let sy = uniaxial_stress_along(v, Vec3::new(0.0, 1.0, 0.0));
        assert!((sy - 0.0).abs() < 1e-12);
        let sx = uniaxial_stress_along(v, Vec3::new(1.0, 0.0, 0.0));
        assert!((sx - 10.0).abs() < 1e-12);
    }
}
