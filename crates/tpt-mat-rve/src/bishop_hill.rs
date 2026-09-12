//! Bishop-Hill (1951) maximum-work / Taylor-factor solver.
//!
//! The Taylor (1938) assumption in polycrystal plasticity is that
//! each grain undergoes the same macroscopic strain `eps`, and the
//! slip systems that minimise the dissipated slip (equivalently
//! maximise the macroscopic work) are activated.  Bishop & Hill
//! (1951) showed that this reduces to a pair of dual linear programs:
//!
//! `min_γ  Σ |γ^α|   s.t.   Σ γ^α P^α = eps`   (primal),
//!
//! `max_σ  σ : eps   s.t.   |σ : P^α| <= τ_c`  (dual),
//!
//! where `P^α` is the Schmid tensor of slip system α.  For FCC,
//! the yield-surface vertices are the intersections of 5 or 6 of
//! the 24 yield facets (the 56 Bishop–Hill corners arise through
//! the sign freedom of the resolved shear on each system).
//!
//! # Implementation
//!
//! The solver enumerates the vertices of the stress yield polytope
//! `|σ : P^α| ≤ τ_c` — five tight slip systems together with a sign
//! pattern — keeps the maximum-work dual-feasible vertex, and recovers
//! the slip rates from the tight set by an exact linear solve.  The
//! primal–dual pair satisfies `σ : ε = τ_c Σ|γ^α|` by LP duality, and
//! `M = Σ|γ^α| / ε_eq`.  For a random FCC polycrystal this recovers
//! `M = 3.06` (Taylor, 1938).
//!
//! [`bishop_hill_lemke`] solves the same LCP through this exact core.
//!
//! # References
//!
//! Reference: Bishop, J. F. W. & Hill, R. (1951).  "A theory of
//! the plastic distortion of a polycrystalline aggregate of pure
//! metals." Phil. Mag. 42, 1298-1307.

use std::ops::Sub;

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

/// Extract uniaxial stress along the unit axis `t` from a Voigt
/// stress `sigma` via `sigma_yy = sigma : (t outer t)`.
pub fn uniaxial_stress_along(sigma: Vec6, t: Vec3) -> f64 {
    let ts = t.sym_outer(t);
    let v = Vec6::from_sym_mat3(ts);
    sigma.double_dot(v)
}

/// Compute the Taylor factor of a *single* grain with the crystal
/// frame aligned such that the tensile axis lies along `tensile_axis`.
///
/// The applied strain is the unit *uniaxial deviator* `(3/2)(t ⊗ t - I/3)`,
/// whose axial component and von-Mises magnitude both equal 1, so the
/// returned `M = Σ |γ^α|` directly mirrors the classical
/// Bishop–Hill values (`[001] → 2.449`, `[111] → 3.674`).
///
/// `slips` is an optional explicit slip-system list (defaults to the
/// FCC `{111}<110>` family); every system shares the same `τ_c`.
pub fn bishop_hill_taylor_factor_axis(
    tensile_axis: [f64; 3],
    tau_c: f64,
    slips: Option<&[SlipSystem]>,
) -> BishopHillResult {
    let t = Vec3::from(tensile_axis).normalized();
    let tt = Vec6::from_sym_mat3(t.sym_outer(t));
    let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    let eps = tt.sub(id.scale(1.0 / 3.0)).scale(1.5);
    bishop_hill_for_strain(eps, tau_c, slips)
}

/// Solve the exact L1 Bishop-Hill problem for a *given* macroscopic
/// strain `eps`.
pub fn bishop_hill_taylor_factor(eps: Vec6, tau_c: f64) -> BishopHillResult {
    bishop_hill_for_strain(eps, tau_c, None)
}

/// Solve the Bishop–Hill maximum-work problem for the FCC
/// `{111}<110>` slip family — the LCP `min Σ|γ^α|` s.t.
/// `Σ γ^α P^α = ε'`, `γ^α ≥ 0` — via exact primal–dual vertex
/// enumeration.  For a random FCC polycrystal this recovers the
/// classical `M = 3.06`.
pub fn bishop_hill_lemke(eps: Vec6, tau_c: f64) -> BishopHillResult {
    let slips = CrystalStructure::FCC.slip_systems();
    let schmid = build_schmid(&slips);
    bishop_hill_lemke_with_slips(eps, tau_c, &slips, &schmid)
}

/// Convert 6-component Voigt deviatoric to 5-component [e11, e22, e12, e13, e23]
fn voigt6_to_dev5(v: Vec6) -> [f64; 5] {
    let tr = (v.data[0] + v.data[1] + v.data[2]) / 3.0;
    [
        v.data[0] - tr,
        v.data[1] - tr,
        v.data[5],
        v.data[4],
        v.data[3],
    ]
}

/// Convert 5-component deviatoric back to 6-component Voigt (traceless)
fn dev5_to_voigt6(dev5: &[f64; 5]) -> Vec6 {
    let e11 = dev5[0];
    let e22 = dev5[1];
    let e12 = dev5[2];
    let e13 = dev5[3];
    let e23 = dev5[4];
    let e33 = -e11 - e22;
    Vec6::new(e11, e22, e33, e23, e13, e12)
}

/// Exact LCP core shared by the public entry points: dispatches to the
/// primal–dual vertex enumeration in [`bishop_hill_from_schmid`].
fn bishop_hill_lemke_solve(
    eps: Vec6,
    tau_c: f64,
    slips: &[SlipSystem],
    schmid: &[Vec6],
) -> BishopHillResult {
    let _ = slips;
    bishop_hill_from_schmid(eps, tau_c, schmid.len(), schmid)
}

/// Invert a 5×5 matrix (Gauss-Jordan)
fn invert_5x5(m: &mut [[f64; 5]; 5]) -> Option<[[f64; 5]; 5]> {
    let mut a = *m;
    let mut inv = [[0.0_f64; 5]; 5];
    for i in 0..5 {
        inv[i][i] = 1.0;
    }
    for i in 0..5 {
        // Find pivot
        let mut piv = i;
        for k in (i + 1)..5 {
            if a[k][i].abs() > a[piv][i].abs() {
                piv = k;
            }
        }
        if a[piv][i].abs() < 1e-12 {
            return None;
        }
        if piv != i {
            a.swap(i, piv);
            inv.swap(i, piv);
        }
        let piv_val = a[i][i];
        for j in 0..5 {
            a[i][j] /= piv_val;
            inv[i][j] /= piv_val;
        }
        for k in 0..5 {
            if k != i {
                let factor = a[k][i];
                for j in 0..5 {
                    a[k][j] -= factor * a[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }
    Some(inv)
}

/// Like [`bishop_hill_lemke`] but with caller-supplied slip systems.
pub fn bishop_hill_lemke_with_slips(
    eps: Vec6,
    tau_c: f64,
    slips: &[SlipSystem],
    schmid: &[Vec6],
) -> BishopHillResult {
    bishop_hill_lemke_solve(eps, tau_c, slips, schmid)
}

/// Shared driver: build the Schmid tensors and hand off to the solver.
fn bishop_hill_for_strain(eps: Vec6, tau_c: f64, slips: Option<&[SlipSystem]>) -> BishopHillResult {
    let slips = slips
        .map(<[SlipSystem]>::to_vec)
        .unwrap_or_else(|| CrystalStructure::FCC.slip_systems());
    let schmid = build_schmid(&slips);
    bishop_hill_from_schmid(eps, tau_c, slips.len(), &schmid)
}

/// Upper bound on the number of enumerated candidate active sets
/// before we fall back to the regularised least-squares estimate.
/// `(12 choose 5) = 792` for FCC stays fast; `(48 choose 5) ≈ 1.7 × 10⁶`
/// does not.
const MAX_COMBOS: usize = 200_000;

/// Apply the metric `G` of the 5-component traceless-tensor
/// representation `[11, 22, 12, 13, 23]`, defined by
/// `x : y = x · (G y)` so that the result matches the Voigt double-dot
/// product of the full symmetric tensors.
fn metric_apply(x: &[f64; 5]) -> [f64; 5] {
    [
        2.0 * x[0] + x[1],
        x[0] + 2.0 * x[1],
        2.0 * x[2],
        2.0 * x[3],
        2.0 * x[4],
    ]
}

/// Tensor inner product of two traceless tensors in the 5-component
/// representation (`x : y = x · G y`).
fn dot5(x: &[f64; 5], y: &[f64; 5]) -> f64 {
    let gy = metric_apply(y);
    x.iter().zip(gy).map(|(a, b)| a * b).sum()
}

/// Exact Bishop–Hill solver: minimise `Σ|γ^α|` subject to
/// `Σ γ^α P^α = ε'` over the given traceless Schmid tensors.
///
/// The optimum is found by enumerating the vertices of the dual stress
/// polytope `|σ : P^α| ≤ τ_c` (5 tight systems × 2⁵ sign patterns),
/// keeping the maximum-work vertex, and recovering the slip rates from
/// the tight set.  LP duality guarantees `σ : ε = τ_c Σ|γ^α|` at the
/// optimum, so `M = Σ|γ^α| / ε_eq` matches the classical values.
#[allow(clippy::too_many_lines)]
fn bishop_hill_from_schmid(eps: Vec6, tau_c: f64, n: usize, schmid: &[Vec6]) -> BishopHillResult {
    let eps_dev = deviator(eps);
    let eps_eq = strain_von_mises(eps);
    if n < 5 || eps_eq < 1e-15 || tau_c.abs() < 1e-15 || binom(n, 5) > MAX_COMBOS {
        return bishop_hill_lsq_fallback(eps, tau_c, schmid);
    }

    // Deviatoric Schmid tensors and applied strain in the 5-component
    // representation.
    let p5: Vec<[f64; 5]> = schmid.iter().map(|&p| voigt6_to_dev5(p)).collect();
    let e5 = voigt6_to_dev5(eps_dev);
    let tc = tau_c.abs();
    let feasibility_tol = tc * (1.0 + 1.0e-7);

    // Pass 1 — maximum work over all dual-feasible stress vertices.
    let mut best_work = f64::NEG_INFINITY;
    for combo in combinations(n, 5) {
        // Row `i` is the metric-weighted Schmid vector of `combo[i]`, so
        // that `σ : P^α = A[i] · σ5` for the 5-component stress `σ5`.
        let mut a = [[0.0_f64; 5]; 5];
        for (i, &sys) in combo.iter().enumerate() {
            a[i] = metric_apply(&p5[sys]);
        }
        // Sign patterns only flip right-hand-side signs, so a single
        // inversion per subset covers all 2⁵ patterns.
        let Some(inv) = invert_5x5(&mut a) else {
            continue;
        };
        for pattern in 0..32_usize {
            let mut sigma5 = [0.0_f64; 5];
            for (i, s) in sigma5.iter_mut().enumerate() {
                let mut acc = 0.0_f64;
                for (j, item) in inv[i].iter().enumerate() {
                    let sign = if (pattern >> j) & 1 == 0 { tc } else { -tc };
                    acc += item * sign;
                }
                *s = acc;
            }
            // Dual feasibility: `|σ : P^β| ≤ τ_c` for every system.
            if p5.iter().any(|p| dot5(&sigma5, p).abs() > feasibility_tol) {
                continue;
            }
            let work = dot5(&sigma5, &e5);
            if work > best_work {
                best_work = work;
            }
        }
    }
    if !best_work.is_finite() {
        // No feasible vertex: strain not reconstructible by this slip set.
        return bishop_hill_lsq_fallback(eps, tau_c, schmid);
    }

    // Pass 2 — among the (near-)tied optimal vertices, recover the slip
    // rates from the tight set and keep the most fully determined one.
    let mut best: Option<(usize, BishopHillResult)> = None;
    for combo in combinations(n, 5) {
        let mut a = [[0.0_f64; 5]; 5];
        for (i, &sys) in combo.iter().enumerate() {
            a[i] = metric_apply(&p5[sys]);
        }
        let Some(inv) = invert_5x5(&mut a) else {
            continue;
        };
        for pattern in 0..32_usize {
            let mut sigma5 = [0.0_f64; 5];
            for (i, s) in sigma5.iter_mut().enumerate() {
                let mut acc = 0.0_f64;
                for (j, item) in inv[i].iter().enumerate() {
                    let sign = if (pattern >> j) & 1 == 0 { tc } else { -tc };
                    acc += item * sign;
                }
                *s = acc;
            }
            if p5.iter().any(|p| dot5(&sigma5, p).abs() > feasibility_tol) {
                continue;
            }
            let work = dot5(&sigma5, &e5);
            if work < best_work - 1.0e-9 * tc * eps_eq {
                continue;
            }
            // Tight systems: `|σ : P^β| ≈ τ_c`, with the yield-plane sign.
            let mut tight = Vec::new();
            let mut signs = Vec::new();
            for (sys, p) in p5.iter().enumerate() {
                let d = dot5(&sigma5, p);
                if d.abs() >= tc * (1.0 - 1.0e-6) {
                    tight.push(sys);
                    signs.push(d.signum());
                }
            }
            let Some((rates, n_active)) = recover_slip_rates(&p5, &tight, &signs, &e5) else {
                continue;
            };
            let improves = match &best {
                None => true,
                Some((na, _)) => n_active > *na,
            };
            if !improves {
                continue;
            }
            let mut plastic = [0.0_f64; 6];
            let mut slip_full = vec![0.0_f64; n];
            let mut sum_abs = 0.0_f64;
            for &(sys, g) in &rates {
                slip_full[sys] = g;
                sum_abs += g.abs();
                for (j, pl) in plastic.iter_mut().enumerate() {
                    *pl += g * schmid[sys].data[j];
                }
            }
            // Complementary slackness: a recovered primal lies on the
            // optimal face only if its cost matches the vertex work
            // (`sigma : eps = tau_c * Sum |gamma|`).
            if (sum_abs * tc - work).abs() > 1.0e-6 * work.abs().max(1.0) {
                continue;
            }
            let result = BishopHillResult {
                n_active,
                slip_rates: slip_full,
                plastic_strain: Vec6::new(
                    plastic[0], plastic[1], plastic[2], plastic[3], plastic[4], plastic[5],
                ),
                stress: dev5_to_voigt6(&sigma5),
                taylor_factor: sum_abs / eps_eq,
            };
            best = Some((n_active, result));
        }
    }

    match best {
        Some((_, result)) => result,
        None => bishop_hill_lsq_fallback(eps, tau_c, schmid),
    }
}

/// Recover non-negative slip rates on the tight systems of an optimal
/// stress vertex: find `gamma >= 0` with `Sum_i gamma_i s_i P^(a_i) = eps'`
/// by solving every 5-subset of the tight set exactly.  Returns the
/// signed rates together with the number of strictly active systems.
fn recover_slip_rates(
    p5: &[[f64; 5]],
    tight: &[usize],
    signs: &[f64],
    e5: &[f64; 5],
) -> Option<(Vec<(usize, f64)>, usize)> {
    if tight.len() < 5 || tight.len() > 12 {
        return None;
    }
    let mut best: Option<(usize, Vec<(usize, f64)>)> = None;
    for combo in combinations(tight.len(), 5) {
        // Column `i` is the signed Schmid vector `s_i P^(a_i)`.
        let mut b = [[0.0_f64; 5]; 5];
        for (i, &ci) in combo.iter().enumerate() {
            let col = p5[tight[ci]];
            for (r, item) in b.iter_mut().enumerate() {
                item[i] = signs[ci] * col[r];
            }
        }
        let Some(inv) = invert_5x5(&mut b) else {
            continue;
        };
        let mut g = [0.0_f64; 5];
        for (i, gi) in g.iter_mut().enumerate() {
            *gi = inv[i].iter().zip(e5).map(|(a, &e)| a * e).sum();
        }
        if g.iter().any(|&x| x < -1.0e-9) {
            continue;
        }
        let n_active = g.iter().filter(|&&x| x > 1.0e-9).count();
        let improves = match &best {
            None => true,
            Some((na, _)) => n_active > *na,
        };
        if improves {
            best = Some((
                n_active,
                combo
                    .iter()
                    .copied()
                    .zip(g)
                    .map(|(ci, gi)| (tight[ci], signs[ci] * gi))
                    .collect(),
            ));
        }
    }
    // Degenerate vertices (more than 5 tight systems): the 5-subset
    // solves above carry structural zeros.  Also try the (regularised)
    // minimum-norm solution over the full tight set, `gamma =
    // B^T (B B^T + eps I)^-1 e5`, which recovers the symmetric
    // multi-system Taylor state for `[111]`-type axes.  Candidates whose
    // cost does not match the vertex work are rejected by the caller.
    if tight.len() > 5 {
        let k = tight.len();
        // B: 5 x k with columns `signs[i] * p5[tight[i]]`.
        let mut bbT = [[0.0_f64; 5]; 5];
        let mut b_e5 = [0.0_f64; 5];
        for i in 0..k {
            let col = p5[tight[i]];
            for (r, item) in bbT.iter_mut().enumerate() {
                for c in 0..5 {
                    item[c] += signs[i] * col[r] * signs[i] * col[c];
                }
            }
            for (r, item) in b_e5.iter_mut().enumerate() {
                *item += signs[i] * col[r] * e5[r];
            }
        }
        // Ridge for rank-deficient tight sets (validated by the residual).
        let max_diag = bbT
            .iter()
            .enumerate()
            .map(|(i, row)| row[i].abs())
            .fold(0.0_f64, f64::max);
        if max_diag > 0.0 {
            for (i, row) in bbT.iter_mut().enumerate() {
                row[i] += 1.0e-9 * max_diag;
            }
            if let Some(inv) = invert_5x5(&mut bbT) {
                // x = (B B^T + eps I)^-1 e5, then gamma = B^T x.
                let mut x = [0.0_f64; 5];
                for (i, xi) in x.iter_mut().enumerate() {
                    *xi = inv[i].iter().zip(b_e5).map(|(a, e)| a * e).sum();
                }
                let mut g = vec![0.0_f64; k];
                let mut rec = [0.0_f64; 5];
                for i in 0..k {
                    let col = p5[tight[i]];
                    g[i] = signs[i] * col.iter().zip(x).map(|(a, b)| a * b).sum::<f64>();
                    for (r, item) in rec.iter_mut().enumerate() {
                        *item += signs[i] * g[i] * col[r];
                    }
                }
                let res: f64 = rec.iter().zip(e5).map(|(a, &b)| (a - b) * (a - b)).sum();
                let target: f64 = e5.iter().map(|&x| x * x).sum::<f64>().max(1e-30);
                if g.iter().all(|&x| x > -1.0e-9) && res <= 1.0e-9 * target {
                    let n_active = g.iter().filter(|&&x| x > 1.0e-9).count();
                    let improves = match &best {
                        None => true,
                        Some((na, _)) => n_active > *na,
                    };
                    if improves {
                        best = Some((
                            n_active,
                            tight
                                .iter()
                                .copied()
                                .enumerate()
                                .map(|(i, sys)| (sys, signs[i] * g[i]))
                                .collect(),
                        ));
                    }
                }
            }
        }
    }
    best.map(|(n_active, rates)| (rates, n_active))
}

/// Deviatoric part of a Voigt strain.
fn deviator(v: Vec6) -> Vec6 {
    let tr = (v[0] + v[1] + v[2]) / 3.0;
    Vec6::new(v[0] - tr, v[1] - tr, v[2] - tr, v[3], v[4], v[5])
}

/// All `k`-subsets of `{0..n}` (lexicographic).
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    if k > n {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut combo: Vec<usize> = (0..k).collect();
    loop {
        out.push(combo.clone());
        if let Some(i) = (0..k).rev().find(|&i| combo[i] + (k - i) < n) {
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

/// Binomial coefficient `C(n, k)`.
fn binom(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1usize;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

/// Regularised least-squares slip estimate, kept as a fallback for
/// degenerate or oversized slip sets that cannot be solved exactly.
fn bishop_hill_lsq_fallback(eps: Vec6, tau_c: f64, schmid: &[Vec6]) -> BishopHillResult {
    let n = schmid.len();
    let p_pinv = pseudo_inverse_schmid(schmid);
    let mut gamma = vec![0.0_f64; n];
    for i in 0..n {
        let mut acc = 0.0;
        for j in 0..6 {
            acc += p_pinv[i][j] * eps[j];
        }
        gamma[i] = acc;
    }
    let m = taylor_factor_from_gamma(&gamma, &eps);
    let mut plastic = [0.0_f64; 6];
    for (i, g) in gamma.iter().enumerate() {
        for j in 0..6 {
            plastic[j] += g * schmid[i].data[j];
        }
    }
    let _ = tau_c;
    BishopHillResult {
        n_active: gamma.iter().filter(|&&g| g.abs() > 1e-9).count(),
        slip_rates: gamma,
        plastic_strain: Vec6::new(
            plastic[0], plastic[1], plastic[2], plastic[3], plastic[4], plastic[5],
        ),
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

/// von-Mises equivalent strain of a Voigt strain tensor.
///
/// `eps_eq = sqrt(2/3 ‖ε'‖_F²)` where `ε'` is the deviatoric part;
/// for a unit uniaxial deviator strain this equals the axial strain,
/// matching the classical convention used for the Bishop–Hill factors.
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
    let j2 = 0.5
        * (e_dev[0] * e_dev[0]
            + e_dev[1] * e_dev[1]
            + e_dev[2] * e_dev[2]
            + 2.0 * e_dev[3] * e_dev[3]
            + 2.0 * e_dev[4] * e_dev[4]
            + 2.0 * e_dev[5] * e_dev[5]);
    // j2 = ½‖ε'‖_F² ⇒ eps_eq = sqrt(2/3 · 2 j2) = sqrt(4/3 j2).
    ((4.0 / 3.0) * j2).sqrt()
}

/// Pseudo-inverse of the 6 x n Schmid matrix `P` via
/// `P^+ = P^T (P P^T + eps I)^{-1}` (true Moore–Penrose pseudo-inverse
/// for rank-deficient `P`).
fn pseudo_inverse_schmid(schmid: &[Vec6]) -> Vec<Vec<f64>> {
    let n = schmid.len();
    let mut ppt = vec![vec![0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..n {
                ppt[i][j] += schmid[k].data[i] * schmid[k].data[j];
            }
        }
    }
    let diag = (0..6).map(|i| ppt[i][i]).fold(0.0_f64, f64::max);
    let eps_reg = diag * 1e-6;
    for i in 0..6 {
        ppt[i][i] += eps_reg;
    }
    let ppt_inv = invert_kxk_pub(&mut ppt).unwrap_or_else(|| vec![vec![0.0; 6]; 6]);
    let mut pinv = vec![vec![0.0_f64; 6]; n];
    for i in 0..n {
        for j in 0..6 {
            let mut acc = 0.0;
            for k in 0..6 {
                acc += schmid[i].data[k] * ppt_inv[k][j];
            }
            pinv[i][j] = acc;
        }
    }
    pinv
}

fn build_schmid(slips: &[SlipSystem]) -> Vec<Vec6> {
    slips
        .iter()
        .map(|s| {
            let p = SlipSystem::schmid_tensor(s.slip_direction, s.plane_normal);
            let mut v = Vec6::from_sym_mat3(p);
            // Ensure traceless (Schmid tensor should be traceless: s·n = 0)
            let tr = (v.data[0] + v.data[1] + v.data[2]) / 3.0;
            v.data[0] -= tr;
            v.data[1] -= tr;
            v.data[2] -= tr;
            v
        })
        .collect()
}

fn solve_slip_rates(schmid: &[Vec6], combo: &[usize], eps: Vec6) -> Vec<f64> {
    let k = combo.len();
    if k == 0 {
        return vec![];
    }
    if k == 5 {
        // Exact 5x5 system for the 5 independent deviatoric strain components
        // Using factors of 2 for shear components (Voigt convention)
        let mut m = vec![vec![0.0_f64; 5]; 5];
        let mut b = vec![0.0_f64; 5];
        for (i, &a) in combo.iter().enumerate() {
            let p = &schmid[a].data;
            m[i][0] = p[0] - p[2];
            m[i][1] = p[1] - p[2];
            m[i][2] = 2.0 * p[3];
            m[i][3] = 2.0 * p[4];
            m[i][4] = 2.0 * p[5];
        }
        b[0] = eps.data[0] - eps.data[2];
        b[1] = eps.data[1] - eps.data[2];
        b[2] = 2.0 * eps.data[3];
        b[3] = 2.0 * eps.data[4];
        b[4] = 2.0 * eps.data[5];
        match invert_kxk(&mut m) {
            Some(inv) => {
                let mut g = vec![0.0_f64; 5];
                for i in 0..5 {
                    let mut acc = 0.0;
                    for j in 0..5 {
                        acc += inv[i][j] * b[j];
                    }
                    g[i] = acc;
                }
                g
            }
            None => vec![0.0; 5],
        }
    } else {
        let mut ptp = vec![vec![0.0_f64; k]; k];
        let mut pt_eps = vec![0.0_f64; k];
        for (i, &a) in combo.iter().enumerate() {
            for (j, &b_idx) in combo.iter().enumerate() {
                ptp[i][j] = schmid[a].double_dot(schmid[b_idx]);
            }
            pt_eps[i] = schmid[a].double_dot(eps);
        }
        let mut ptp_reg = ptp.clone();
        let diag = (0..k).map(|i| ptp[i][i]).fold(0.0_f64, f64::max);
        let eps_reg = diag * 1e-6;
        for i in 0..k {
            ptp_reg[i][i] += eps_reg;
        }
        let result = invert_kxk_pub(&mut ptp_reg);
        if let Some(inv) = result {
            let mut g = vec![0.0_f64; k];
            for i in 0..k {
                let mut acc = 0.0;
                for j in 0..k {
                    acc += inv[i][j] * pt_eps[j];
                }
                g[i] = acc;
            }
            g
        } else if let Some(inv) = invert_kxk_pub(&mut ptp) {
            let mut g = vec![0.0_f64; k];
            for i in 0..k {
                let mut acc = 0.0;
                for j in 0..k {
                    acc += inv[i][j] * pt_eps[j];
                }
                g[i] = acc;
            }
            g
        } else {
            vec![0.0; k]
        }
    }
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
            eprintln!(
                "DEBUG invert_kxk: pivot[{i}] = {} < 1e-6, matrix singular",
                a[piv][i]
            );
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn debug_elastic_110_min_residual() {
        let slips = CrystalStructure::FCC.slip_systems();
        let schmid = build_schmid(&slips);
        let t = Vec3::new(1.0, 1.0, 0.0).normalized();
        let tt = Vec6::from_sym_mat3(t.sym_outer(t));
        let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        // Elastic strain for uniaxial tension along [110] with ν=0.5:
        // ε = (1/E)(σ - 0.5 tr(σ) I) = (σ_0/E)(t⊗t - 0.5 I)
        let eps = tt.sub(id.scale(0.5)).scale(1.0);
        let eps_dev = deviator(eps);
        let residual_scale = eps_dev.double_dot(eps_dev);

        let mut min5 = f64::INFINITY;
        let mut min_combo = None;
        for combo in combinations(12, 5) {
            let gamma = solve_slip_rates(&schmid, &combo, eps_dev);
            if gamma.iter().all(|&g| g.abs() < 1e-15) {
                continue;
            }
            let mut rec = [0.0_f64; 6];
            for (i, &a) in combo.iter().enumerate() {
                for j in 0..6 {
                    rec[j] += gamma[i] * schmid[a].data[j];
                }
            }
            let diff = Vec6::new(rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]).sub(eps_dev);
            let r = diff.double_dot(diff);
            if r < min5 {
                min5 = r;
                min_combo = Some(combo.clone());
            }
        }
        eprintln!(
            "[110] elastic min 5-residual={} (tol={})",
            min5,
            1e-4 * residual_scale
        );
        if let Some(combo) = min_combo {
            eprintln!("[110] elastic min 5-set combo: {:?}", combo);
        }
    }

    #[test]
    fn debug_110_six_combo_residual() {
        let slips = CrystalStructure::FCC.slip_systems();
        let schmid = build_schmid(&slips);
        let t = Vec3::new(1.0, 1.0, 0.0).normalized();
        let tt = Vec6::from_sym_mat3(t.sym_outer(t));
        let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        let eps = tt.sub(id.scale(1.0 / 3.0)).scale(1.5);
        let eps_dev = deviator(eps);

        let mut min6 = f64::INFINITY;
        let mut min6_combo = None;
        for combo in combinations(12, 6) {
            let gamma = solve_slip_rates(&schmid, &combo, eps_dev);
            if gamma.iter().all(|&g| g.abs() < 1e-15) {
                continue;
            }
            // All active slip rates should have the same sign (for a valid vertex)
            let first_sign = gamma[0].signum();
            if gamma.iter().any(|&g| g.signum() != first_sign) {
                continue;
            }
            let mut rec = [0.0_f64; 6];
            for (i, &a) in combo.iter().enumerate() {
                for j in 0..6 {
                    rec[j] += gamma[i] * schmid[a].data[j];
                }
            }
            let diff = Vec6::new(rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]).sub(eps_dev);
            let r = diff.double_dot(diff);
            if r < min6 {
                min6 = r;
                min6_combo = Some(combo.clone());
            }
        }
        eprintln!(
            "[110] unit min 6-residual={} (tol={})",
            min6,
            1e-4 * eps_dev.double_dot(eps_dev)
        );
        if let Some(combo) = min6_combo {
            eprintln!("[110] unit min 6-set combo: {:?}", combo);
        }
    }

    #[test]
    fn debug_110_exact_solver() {
        let tau_c = 1.0;
        let r = bishop_hill_taylor_factor_axis([1.0, 1.0, 0.0], tau_c, None);
        eprintln!("[110] M = {}", r.taylor_factor);
        eprintln!("[110] stress = {:?}", r.stress.data);
        eprintln!("[110] n_active = {}", r.n_active);
        eprintln!("[110] gamma = {:?}", r.slip_rates);
    }

    #[test]
    fn debug_random_min_residual() {
        let slips = CrystalStructure::FCC.slip_systems();
        let schmid = build_schmid(&slips);
        let mut rng = 0xDEADBEEFu64;
        let mut next = || {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            rng as f64 / u64::MAX as f64
        };
        let mut worst_residual = 0.0;
        let mut worst_axis = None;
        for _ in 0..100 {
            let axis = loop {
                let u = 2.0 * next() - 1.0;
                let v = 2.0 * next() - 1.0;
                let s = u * u + v * v;
                if s < 1.0 && s > 0.0 {
                    let factor = 2.0 * (1.0 - s).sqrt();
                    break [u * factor, v * factor, 1.0 - 2.0 * s];
                }
            };
            let t = Vec3::from(axis).normalized();
            let tt = Vec6::from_sym_mat3(t.sym_outer(t));
            let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
            let eps = tt.sub(id.scale(1.0 / 3.0)).scale(1.5);
            let eps_dev = deviator(eps);
            let mut min_residual = f64::INFINITY;
            let mut count_small = 0;
            for combo in combinations(12, 5) {
                let gamma = solve_slip_rates(&schmid, &combo, eps_dev);
                if gamma.iter().all(|&g| g.abs() < 1e-15) {
                    continue;
                }
                let mut rec = [0.0_f64; 6];
                for (i, &a) in combo.iter().enumerate() {
                    for j in 0..6 {
                        rec[j] += gamma[i] * schmid[a].data[j];
                    }
                }
                let diff = Vec6::new(rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]).sub(eps_dev);
                let r = diff.double_dot(diff);
                if r < min_residual {
                    min_residual = r;
                }
                if r < 1e-1 * eps_dev.double_dot(eps_dev) {
                    count_small += 1;
                }
            }
            if min_residual > worst_residual {
                worst_residual = min_residual;
                worst_axis = Some((axis, count_small));
            }
        }
        eprintln!(
            "worst min residual over 100 random axes: {} (tol = {})",
            worst_residual,
            1e-4 * 1.5
        );
        if let Some((axis, count)) = worst_axis {
            eprintln!("worst axis: {:?}, count_small@0.1: {}", axis, count);
        }
    }

    #[test]
    fn fcc_benchmark_axes_recover_classical_taylor_factors() {
        let tau_c = 1.0;
        let m100 = bishop_hill_taylor_factor_axis([1.0, 0.0, 0.0], tau_c, None);
        let m110 = bishop_hill_taylor_factor_axis([1.0, 1.0, 0.0], tau_c, None);
        let sqrt6 = 6.0_f64.sqrt();
        let three_sqrt6_over2 = 1.5 * sqrt6;
        assert!(
            (m100.taylor_factor - sqrt6).abs() < 1e-3,
            "[001] M = {} != √6 = {}",
            m100.taylor_factor,
            sqrt6
        );
        assert!(
            (m110.taylor_factor - three_sqrt6_over2).abs() < 1e-3,
            "[110] M = {} != 3√6/2 = {}",
            m110.taylor_factor,
            three_sqrt6_over2
        );
        // [111] needs a 6-active-system vertex; the vertex enumeration
        // covers it through its 5-subsets (degenerate vertex).
        let m111 = bishop_hill_taylor_factor_axis([1.0, 1.0, 1.0], tau_c, None);
        assert!(
            (m111.taylor_factor - three_sqrt6_over2).abs() < 1e-3,
            "[111] M = {} != 3√6/2 = {}",
            m111.taylor_factor,
            three_sqrt6_over2
        );
        // All three vertices must be feasible: |σ : P^β| ≤ τ_c.
        for r in [&m100, &m110, &m111] {
            assert!(r.stress != Vec6::ZERO);
            // The Taylor stress is multi-critical: several slip systems
            // sit on the yield surface simultaneously (5-8 for these
            // axes; the exact count depends on the vertex chosen).
            let n_at_yield = CrystalStructure::FCC
                .slip_systems()
                .iter()
                .map(|s| {
                    let p = SlipSystem::schmid_tensor(s.slip_direction, s.plane_normal);
                    r.stress.double_dot(Vec6::from_sym_mat3(p)).abs()
                })
                .filter(|&tau| (tau - tau_c).abs() < 1.0e-6)
                .count();
            assert!(
                n_at_yield >= 5,
                "only {n_at_yield} systems at yield (expected >= 5)"
            );
            assert!(
                r.n_active >= 3,
                "recovered slip state has only {} active systems",
                r.n_active
            );
        }
    }

    #[test]
    fn fcc_random_average_taylor_factor_is_about_3_06() {
        let tau_c = 1.0;
        let mut rng = 0xDEADBEEFu64;
        let mut next = || {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            rng as f64 / u64::MAX as f64
        };
        let mut sum = 0.0;
        let n_samples = 1024;
        for i in 0..n_samples {
            let axis = loop {
                let u = 2.0 * next() - 1.0;
                let v = 2.0 * next() - 1.0;
                let s = u * u + v * v;
                if s < 1.0 && s > 0.0 {
                    let factor = 2.0 * (1.0 - s).sqrt();
                    break [u * factor, v * factor, 1.0 - 2.0 * s];
                }
            };
            let t = Vec3::from(axis).normalized();
            let tt = Vec6::from_sym_mat3(t.sym_outer(t));
            let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
            let eps = tt.sub(id.scale(1.0 / 3.0)).scale(1.5);
            let r = bishop_hill_lemke(eps, tau_c);
            if !r.taylor_factor.is_finite() || r.taylor_factor <= 0.0 {
                eprintln!(
                    "FAIL at sample {}: axis={:?}, M={}",
                    i, axis, r.taylor_factor
                );
            }
            sum += r.taylor_factor;
        }
        let mean = sum / n_samples as f64;
        assert!(
            (mean - 3.06).abs() < 0.05,
            "random FCC average M = {mean} != 3.06 ± 0.05"
        );
    }

    #[test]
    fn primal_dual_consistency_max_work_equals_tau_sum_abs() {
        // For the returned vertex: σ : ε = τ_c Σ|γ^α| (LP duality).
        let tau_c = 1.0;
        let r = bishop_hill_taylor_factor_axis([1.0, 2.0, 3.0], tau_c, None);
        assert!(r.stress != Vec6::ZERO);
        let work = r.stress.double_dot(r.plastic_strain);
        let sum_abs: f64 = r.slip_rates.iter().map(|g| g.abs()).sum();
        let rel = (work - tau_c * sum_abs).abs() / (tau_c * sum_abs).max(1e-30);
        assert!(rel < 1e-9, "σ:ε = {work}, τ_c Σ|γ| = {}", tau_c * sum_abs);
    }

    #[test]
    fn plastic_strain_reproduces_deviatoric_input() {
        let axis = [0.0, 0.0, 1.0];
        let t = Vec3::from(axis).normalized();
        let tt = Vec6::from_sym_mat3(t.sym_outer(t));
        let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        let eps = tt.sub(id.scale(1.0 / 3.0)).scale(1.5);
        let r = bishop_hill_lemke(eps, 1.0);
        for j in 0..6 {
            assert!(
                (r.plastic_strain[j] - eps[j]).abs() < 1e-9,
                "plastic[{j}] = {} != eps = {}",
                r.plastic_strain[j],
                eps[j]
            );
        }
    }

    #[test]
    fn stress_is_tensile_along_axis() {
        let tau_c = 1.0;
        let r = bishop_hill_taylor_factor_axis([0.0, 0.0, 1.0], tau_c, None);
        let sy = uniaxial_stress_along(r.stress, Vec3::new(0.0, 0.0, 1.0));
        assert!(sy > 0.0, "uniaxial stress must be tensile, got {sy}");
        let work = r.stress.double_dot(r.plastic_strain);
        let sum_abs: f64 = r.slip_rates.iter().map(|g| g.abs()).sum();
        let m_from_work = work.abs() / (tau_c * sum_abs.max(1e-30));
        assert!(
            (m_from_work - 1.0).abs() < 1e-6,
            "LP duality violated: M_from_work = {}",
            m_from_work
        );
    }

    #[test]
    fn bishop_hill_fcc_random_lemke_returns_finite_result() {
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
            let tt = Vec6::from_sym_mat3(tt_sym);
            let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
            let eps = tt.sub(id.scale(1.0 / 3.0)).scale(1.5);
            let r = bishop_hill_lemke(eps, tau_c);
            assert!(r.taylor_factor.is_finite());
            assert!(r.taylor_factor > 0.0);
        }
    }

    #[test]
    fn bishop_hill_fcc_single_axis_tension_returns_finite_m() {
        let tau_c = 1.0;
        let r = bishop_hill_taylor_factor_axis([0.0, 0.0, 1.0], tau_c, None);
        assert!(r.taylor_factor.is_finite());
        assert!((r.taylor_factor - 6.0_f64.sqrt()).abs() < 1e-2);
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
        assert_eq!(combinations(5, 5).len(), 1);
        assert_eq!(combinations(12, 5).len(), 792);
        assert_eq!(combinations(4, 2).len(), 6);
        assert_eq!(binom(12, 5), 792);
    }

    #[test]
    fn debug_110_best_combo() {
        let slips = CrystalStructure::FCC.slip_systems();
        let schmid = build_schmid(&slips);
        let t = Vec3::new(1.0, 1.0, 0.0).normalized();
        let tt = Vec6::from_sym_mat3(t.sym_outer(t));
        let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        let eps = tt.sub(id.scale(1.0 / 3.0)).scale(1.5);
        let eps_dev = deviator(eps);
        let combo = vec![0usize, 3, 4, 10, 11];
        let gamma = solve_slip_rates(&schmid, &combo, eps_dev);
        eprintln!("[110] combo={:?}", combo);
        eprintln!("[110] gamma={:?}", gamma);
        let mut rec = [0.0_f64; 6];
        for (i, &a) in combo.iter().enumerate() {
            eprintln!("[110] P[{}]={:?}", a, schmid[a].data);
            for j in 0..6 {
                rec[j] += gamma[i] * schmid[a].data[j];
            }
        }
        let diff = Vec6::new(rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]).sub(eps_dev);
        eprintln!("[110] rec={:?}", rec);
        eprintln!("[110] diff={:?}", diff.data);
        eprintln!("[110] residual={}", diff.double_dot(diff));

        let mut ptp = vec![vec![0.0_f64; 5]; 5];
        for (i, &a) in combo.iter().enumerate() {
            for (j, &b) in combo.iter().enumerate() {
                ptp[i][j] = schmid[a].double_dot(schmid[b]);
            }
        }
        eprintln!("[110] PTP det={}", ptp_det(&ptp));
        eprintln!("[110] PTP={:?}", ptp);
    }

    #[test]
    fn test_5x5_invert_matches_identity() {
        let mut m = vec![
            vec![2.0, -1.0, 0.0, 0.0, 0.0],
            vec![-1.0, 2.0, -1.0, 0.0, 0.0],
            vec![0.0, -1.0, 2.0, -1.0, 0.0],
            vec![0.0, 0.0, -1.0, 2.0, -1.0],
            vec![0.0, 0.0, 0.0, -1.0, 2.0],
        ];
        let inv = invert_kxk(&mut m).expect("invert");
        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut g = vec![0.0; 5];
        for i in 0..5 {
            let mut acc = 0.0;
            for j in 0..5 {
                acc += inv[i][j] * b[j];
            }
            g[i] = acc;
        }
        let expected = vec![
            5.833333333333331,
            10.666666666666666,
            13.5,
            13.333333333333332,
            9.166666666666666,
        ];
        for i in 0..5 {
            assert!(
                (g[i] - expected[i]).abs() < 1e-6,
                "g[{i}] = {} != {}",
                g[i],
                expected[i]
            );
        }
    }

    #[test]
    fn test_2x2_invert() {
        let mut m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let inv = invert_kxk(&mut m).expect("invert");
        let expected = vec![vec![-2.0, 1.0], vec![1.5, -0.5]];
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (inv[i][j] - expected[i][j]).abs() < 1e-9,
                    "inv[{i}][{j}] = {} != {}",
                    inv[i][j],
                    expected[i][j]
                );
            }
        }
    }

    #[test]
    fn test_5x5_schmid_system_inverse() {
        let slips = CrystalStructure::FCC.slip_systems();
        let schmid = build_schmid(&slips);
        let combo = vec![0usize, 3, 4, 10, 11];
        let mut m = vec![vec![0.0_f64; 5]; 5];
        for (i, &a) in combo.iter().enumerate() {
            let p = &schmid[a].data;
            m[i][0] = p[0] - p[2];
            m[i][1] = p[1] - p[2];
            m[i][2] = 2.0 * p[3];
            m[i][3] = 2.0 * p[4];
            m[i][4] = 2.0 * p[5];
        }
        let m_orig = m.clone();
        let inv = invert_kxk(&mut m).expect("invert");
        let mut prod = vec![vec![0.0_f64; 5]; 5];
        for i in 0..5 {
            for j in 0..5 {
                let mut acc = 0.0;
                for k in 0..5 {
                    acc += m_orig[i][k] * inv[k][j];
                }
                prod[i][j] = acc;
            }
        }
        for i in 0..5 {
            for j in 0..5 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (prod[i][j] - expected).abs() < 1e-6,
                    "prod[{i}][{j}] = {} != {}",
                    prod[i][j],
                    expected
                );
            }
        }
    }

    #[test]
    fn test_5x5_schmid_system_solve() {
        let slips = CrystalStructure::FCC.slip_systems();
        let schmid = build_schmid(&slips);
        let combo = vec![0usize, 3, 4, 10, 11];
        let t = Vec3::new(1.0, 1.0, 0.0).normalized();
        let tt = Vec6::from_sym_mat3(t.sym_outer(t));
        let id = Vec6::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        let eps = tt.sub(id.scale(1.0 / 3.0)).scale(1.5);
        let eps_dev = deviator(eps);

        let mut m = vec![vec![0.0_f64; 5]; 5];
        let mut b = vec![0.0_f64; 5];
        for (i, &a) in combo.iter().enumerate() {
            let p = &schmid[a].data;
            m[i][0] = p[0] - p[2];
            m[i][1] = p[1] - p[2];
            m[i][2] = 2.0 * p[3];
            m[i][3] = 2.0 * p[4];
            m[i][4] = 2.0 * p[5];
        }
        b[0] = eps_dev.data[0] - eps_dev.data[2];
        b[1] = eps_dev.data[1] - eps_dev.data[2];
        b[2] = 2.0 * eps_dev.data[3];
        b[3] = 2.0 * eps_dev.data[4];
        b[4] = 2.0 * eps_dev.data[5];

        let m_orig = m.clone();
        let inv = invert_kxk(&mut m).expect("invert");
        let mut g = vec![0.0_f64; 5];
        for i in 0..5 {
            let mut acc = 0.0;
            for j in 0..5 {
                acc += inv[i][j] * b[j];
            }
            g[i] = acc;
        }

        let mut mg = vec![0.0_f64; 5];
        for i in 0..5 {
            let mut acc = 0.0;
            for j in 0..5 {
                acc += m_orig[i][j] * g[j];
            }
            mg[i] = acc;
        }

        eprintln!("gamma={:?}", g);
        eprintln!("m*gamma={:?}", mg);
        eprintln!("b={:?}", b);
        eprintln!(
            "diff={:?}",
            mg.iter()
                .zip(b.iter())
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>()
        );

        for i in 0..5 {
            assert!(
                (mg[i] - b[i]).abs() < 1e-6,
                "m*gamma[{}] = {} != b[{}] = {}",
                i,
                mg[i],
                i,
                b[i]
            );
        }
    }

    fn ptp_det(m: &[Vec<f64>]) -> f64 {
        let n = m.len();
        if n == 1 {
            return m[0][0];
        }
        if n == 2 {
            return m[0][0] * m[1][1] - m[0][1] * m[1][0];
        }
        if n == 3 {
            return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        }
        let mut det = 0.0;
        for j in 0..n {
            let mut sub = vec![vec![0.0_f64; n - 1]; n - 1];
            for i in 1..n {
                let mut col = 0;
                for k in 0..n {
                    if k == j {
                        continue;
                    }
                    sub[i - 1][col] = m[i][k];
                    col += 1;
                }
            }
            det += m[0][j] * ptp_det(&sub) * if j % 2 == 0 { 1.0 } else { -1.0 };
        }
        det
    }
}
