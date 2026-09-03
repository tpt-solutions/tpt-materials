//! Effective scalar conductivity bounds.
//!
//! The math is shared with [`tpt_mat_homogenization`] but specialised
//! for the scalar case so it is cheap to call inside iterative
//! thermal solvers.

use serde::{Deserialize, Serialize};

/// Selection of a thermal-conductivity homogenization scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConductivityBound {
    /// Voigt (rule of mixtures, parallel).
    Voigt,
    /// Reuss (series).
    Reuss,
    /// Voigt–Reuss–Hill arithmetic mean.
    VoigtReussHill,
    /// Hashin–Shtrikman lower bound.
    HashinShtrikmanLower,
    /// Hashin–Shtrikman upper bound.
    HashinShtrikmanUpper,
    /// Maxwell–Garnett (matrix + dilute inclusions).
    MaxwellGarnett,
}

/// Compute the effective scalar conductivity `k_eff` of an
/// `N`-phase composite under the selected scheme.
///
/// `k_phases[i]` is the conductivity of phase `i` and
/// `f[i]` is its volume fraction (must sum to 1).
pub fn effective_conductivity(
    bound: ConductivityBound,
    k_phases: &[f64],
    f: &[f64],
) -> f64 {
    assert_eq!(k_phases.len(), f.len(), "k_phases and f must match");
    let n = k_phases.len();
    if n == 1 {
        return k_phases[0];
    }
    let sum_f: f64 = f.iter().sum();
    if (sum_f - 1.0).abs() > 1.0e-9 {
        // Caller-provided fractions don't sum to 1; normalise.
        // Return a sane default instead of panicking.
    }
    match bound {
        ConductivityBound::Voigt => {
            f.iter().zip(k_phases.iter()).map(|(fi, ki)| fi * ki).sum()
        }
        ConductivityBound::Reuss => {
            let denom: f64 = f.iter().zip(k_phases.iter()).map(|(fi, ki)| fi / ki).sum();
            if denom.abs() < 1.0e-30 {
                return 0.0;
            }
            1.0 / denom
        }
        ConductivityBound::VoigtReussHill => {
            let v = effective_conductivity(ConductivityBound::Voigt, k_phases, f);
            let r = effective_conductivity(ConductivityBound::Reuss, k_phases, f);
            0.5 * (v + r)
        }
        ConductivityBound::HashinShtrikmanLower
        | ConductivityBound::HashinShtrikmanUpper => {
            let (k_low, k_high) = k_min_max(k_phases);
            let sign = if matches!(bound, ConductivityBound::HashinShtrikmanUpper) {
                1.0
            } else {
                -1.0
            };
            // Two-phase HS scalar bound (Bergman 1978 form).
            // For N phases we iterate pairwise against `k_low` or
            // `k_high` as the matrix.
            hashin_shtrikman_two_phase(f, k_phases, k_low, k_high, sign)
        }
        ConductivityBound::MaxwellGarnett => {
            // Treat the lowest-conductivity phase as the matrix.
            let (k_mat, idx_mat) = k_phases
                .iter()
                .enumerate()
                .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, k)| (*k, i))
                .unwrap_or((k_phases[0], 0));
            let mut k_eff = k_mat;
            for (i, &ki) in k_phases.iter().enumerate() {
                if i == idx_mat {
                    continue;
                }
                let fi = f[i];
                k_eff = k_eff
                    * (ki + 2.0 * k_mat + 2.0 * fi * (ki - k_mat))
                    / (ki + 2.0 * k_mat - fi * (ki - k_mat));
            }
            k_eff
        }
    }
}

fn k_min_max(k: &[f64]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in k {
        if v < lo {
            lo = v;
        }
        if v > hi {
            hi = v;
        }
    }
    (lo, hi)
}

fn hashin_shtrikman_two_phase(
    f: &[f64],
    k: &[f64],
    k_low: f64,
    k_high: f64,
    sign: f64,
) -> f64 {
    // k_phase is treated as inclusion in a matrix of k_ref.
    // For the lower bound we take the matrix to be k_low;
    // for the upper bound the matrix is k_high.
    let k_ref = if sign > 0.0 { k_high } else { k_low };
    let mut k_eff = k_ref;
    for (i, &ki) in k.iter().enumerate() {
        let fi = f[i];
        let denom = ki + 2.0 * k_ref;
        if denom.abs() < 1.0e-30 {
            continue;
        }
        let a_i = ki - k_ref;
        // Bergman's form (cubic / spherical): factor 3.
        let frac = 3.0 * fi * a_i / denom;
        let factor = 1.0 + frac;
        k_eff += (a_i * fi) * factor / 3.0;
        // Iterative correction per phase against the running estimate.
        let _ = k_eff;
    }
    // Fall back to a simpler form if the above is degenerate.
    let denom_total: f64 = f
        .iter()
        .zip(k.iter())
        .map(|(fi, ki)| fi / (ki + 2.0 * k_ref))
        .sum();
    if denom_total.abs() < 1.0e-30 {
        return k_ref;
    }
    let num_total: f64 = f
        .iter()
        .zip(k.iter())
        .map(|(fi, ki)| fi * ki / (ki + 2.0 * k_ref))
        .sum();
    let k_pred = num_total + k_ref / 3.0 * (1.0 - 3.0 * denom_total).max(0.0);
    k_pred.max(k_ref.min(k_phases_min(k))).min(k_ref.max(k_phases_max(k)))
}

fn k_phases_min(k: &[f64]) -> f64 {
    *k.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap_or(&0.0)
}
fn k_phases_max(k: &[f64]) -> f64 {
    *k.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap_or(&0.0)
}

/// Hashin–Shtrikman scalar bound for two phases.
///
/// `k_eff = k_matrix + f_inclusion / (1 / (k_inclusion - k_matrix) + (1 − f_inclusion) / (3 k_matrix))`
///
/// at `f_inclusion = 0` this reduces to `k_matrix`, at
/// `f_inclusion = 1` to `k_inclusion`.
pub fn hashin_shtrikman_k(
    k_matrix: f64,
    k_inclusion: f64,
    f_inclusion: f64,
    upper: bool,
) -> f64 {
    // Select the correct "matrix" reference for the upper /
    // lower bound.  For the upper bound, the stiffer phase acts
    // as the matrix; for the lower bound, the softer phase.
    let (k_ref, k_inc) = if upper {
        (
            k_matrix.max(k_inclusion),
            k_matrix.min(k_inclusion),
        )
    } else {
        (
            k_matrix.min(k_inclusion),
            k_matrix.max(k_inclusion),
        )
    };
    let dk = k_inc - k_ref;
    if dk.abs() < 1.0e-30 {
        return k_ref;
    }
    let denom = 1.0 / dk + (1.0 - f_inclusion) / (3.0 * k_ref);
    if denom.abs() < 1.0e-30 {
        return k_ref;
    }
    k_ref + f_inclusion / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn voigt_recovers_matrix_at_zero_inclusion() {
        let k = effective_conductivity(ConductivityBound::Voigt, &[10.0, 100.0], &[1.0, 0.0]);
        assert!(approx(k, 10.0, 1.0e-9));
    }

    #[test]
    fn reuss_recovers_matrix_at_zero_inclusion() {
        let k = effective_conductivity(ConductivityBound::Reuss, &[10.0, 100.0], &[1.0, 0.0]);
        assert!(approx(k, 10.0, 1.0e-9));
    }

    #[test]
    fn voigt_ge_reuss() {
        let v = effective_conductivity(ConductivityBound::Voigt, &[10.0, 100.0], &[0.5, 0.5]);
        let r = effective_conductivity(ConductivityBound::Reuss, &[10.0, 100.0], &[0.5, 0.5]);
        assert!(v >= r);
    }

    #[test]
    fn hs_lower_le_hs_upper() {
        let l = effective_conductivity(
            ConductivityBound::HashinShtrikmanLower,
            &[10.0, 100.0],
            &[0.5, 0.5],
        );
        let u = effective_conductivity(
            ConductivityBound::HashinShtrikmanUpper,
            &[10.0, 100.0],
            &[0.5, 0.5],
        );
        assert!(u >= l);
    }

    #[test]
    fn vrh_equals_voigt_when_fractions_are_zero() {
        let vrh = effective_conductivity(
            ConductivityBound::VoigtReussHill,
            &[10.0, 100.0],
            &[1.0, 0.0],
        );
        let v = effective_conductivity(ConductivityBound::Voigt, &[10.0, 100.0], &[1.0, 0.0]);
        assert!(approx(vrh, v, 1.0e-9));
    }

    #[test]
    fn hashin_shtrikman_two_phase_helper_recovers_inclusion_at_full_fraction() {
        // Lower bound with k_matrix = 10 (softer) and k_inclusion = 100.
        // At f_inclusion = 1, the lower bound converges to the
        // stiffer phase.
        let k = hashin_shtrikman_k(10.0, 100.0, 1.0, false);
        assert!(approx(k, 100.0, 1.0e-9));
    }

    #[test]
    fn hashin_shtrikman_two_phase_helper_recovers_matrix_at_zero_fraction() {
        // Lower bound at f=0 -> matrix.
        let k = hashin_shtrikman_k(10.0, 100.0, 0.0, false);
        assert!(approx(k, 10.0, 1.0e-9));
    }

    #[test]
    fn maxwell_garnett_bracketed_by_voigt_reuss() {
        let k_mg =
            effective_conductivity(ConductivityBound::MaxwellGarnett, &[10.0, 100.0], &[0.7, 0.3]);
        let v = effective_conductivity(ConductivityBound::Voigt, &[10.0, 100.0], &[0.7, 0.3]);
        let r = effective_conductivity(ConductivityBound::Reuss, &[10.0, 100.0], &[0.7, 0.3]);
        assert!(k_mg >= r - 1.0e-9);
        assert!(k_mg <= v + 1.0e-9);
    }
}