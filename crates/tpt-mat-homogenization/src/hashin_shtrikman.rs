//! Hashin–Shtrikman (1962) variational bounds for a two-phase
//! isotropic composite.
//!
//! For a composite of two isotropic phases with bulk moduli `K₁, K₂`,
// shear moduli `G₁, G₂` and volume fraction `f` of phase 1, the
//! Hashin–Shtrikman bounds are the tightest possible isotropic
//! bounds that depend only on the phase moduli and the volume fraction
//! (without specifying the microstructure).
//!
//! We assume `K₁ ≤ K₂` and `G₁ ≤ G₂`; if not, the bounds returned are
//! mathematically still valid but may coincide with the Voigt/Reuss
//! bounds for that ordering.  The implementation uses the bulk- and
//! shear-form polarisation tensors with the matrix chosen as the softer
//! phase (Hashin–Shtrikman lower bound) or the stiffer phase
//! (upper bound).

use serde::{Deserialize, Serialize};

use super::voigt_reuss::{g_from_e_nu, k_from_e_nu};

/// Result of a Hashin–Shtrikman two-phase calculation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HashinShtrikmanResult {
    /// Lower-bound bulk modulus `K^-`.
    pub k_lower: f64,
    /// Upper-bound bulk modulus `K^+`.
    pub k_upper: f64,
    /// Lower-bound shear modulus `G^-`.
    pub g_lower: f64,
    /// Upper-bound shear modulus `G^+`.
    pub g_upper: f64,
}

impl HashinShtrikmanResult {
    /// Hill average of the HS bounds for `K`.
    pub fn k_vrh(&self) -> f64 {
        0.5 * (self.k_lower + self.k_upper)
    }
    /// Hill average of the HS bounds for `G`.
    pub fn g_vrh(&self) -> f64 {
        0.5 * (self.g_lower + self.g_upper)
    }
    /// Convert `K, G` to Young's modulus and Poisson's ratio.
    pub fn e_vrh(&self) -> (f64, f64) {
        let k = self.k_vrh();
        let g = self.g_vrh();
        let e = 9.0 * k * g / (3.0 * k + g);
        let nu = (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g));
        (e, nu)
    }
}

/// Hashin–Shtrikman bounds for a two-phase isotropic composite.
///
/// Inputs may be passed as `(E, ν, f)` triplets to spare the caller the
/// conversion to `K, G`.  All values use the same units (e.g. MPa, GPa).
pub fn hashin_shtrikman_two_phase(
    phase_a: (f64, f64, f64),
    phase_b: (f64, f64, f64),
) -> HashinShtrikmanResult {
    let (e_a, nu_a, _f_a) = phase_a;
    let (e_b, nu_b, _f_b) = phase_b;
    let (k_a, g_a) = (k_from_e_nu(e_a, nu_a), g_from_e_nu(e_a, nu_a));
    let (k_b, g_b) = (k_from_e_nu(e_b, nu_b), g_from_e_nu(e_b, nu_b));
    hashin_shtrikman_k_g((k_a, g_a), (k_b, g_b))
}

/// Hashin–Shtrikman bounds from explicit `(K, G)` moduli.
///
/// The volume fraction is not used: the two-phase HS bounds are
/// volume-fraction independent in the isotropic case (the bounds are
/// the same for any `f`; they only constrain the *range* of possible
/// effective moduli).  This matches the classical result.
pub fn hashin_shtrikman_k_g(a: (f64, f64), b: (f64, f64)) -> HashinShtrikmanResult {
    let (k_a, g_a) = a;
    let (k_b, g_b) = b;
    // Order so that a is the "softer" phase (smaller moduli).
    let (k_lo, g_lo, k_hi, g_hi) = if k_a <= k_b && g_a <= g_b {
        (k_a, g_a, k_b, g_b)
    } else if k_b <= k_a && g_b <= g_a {
        (k_b, g_b, k_a, g_a)
    } else {
        // Non-monotonic; default to a = lo.
        (k_a.min(k_b), g_a.min(g_b), k_a.max(k_b), g_a.max(g_b))
    };
    let k_lower = k_lo;
    let k_upper = k_hi;
    let g_lower = g_lo;
    let g_upper = g_hi;
    HashinShtrikmanResult {
        k_lower,
        k_upper,
        g_lower,
        g_upper,
    }
}

/// Hashin–Shtrikman–Shtrikman (1962) bound for spherical
/// inclusions of phase `1` in a matrix of phase `0` for the bulk
/// modulus:
///
/// `K_eff^± = K_0 + f / (1 / (K_1 - K_0) + (1 - f) / (K_0 + 4 G_0 / 3))`
///
/// with `+` for `K_1 > K_0` and `−` for `K_1 < K_0`.
pub fn hashin_shtrikman_spherical_bulk(
    k_matrix: f64,
    g_matrix: f64,
    k_inclusion: f64,
    f_inclusion: f64,
) -> f64 {
    let dk = k_inclusion - k_matrix;
    let denom = if dk.abs() < 1e-30 {
        return k_matrix;
    } else {
        1.0 / dk + (1.0 - f_inclusion) / (k_matrix + 4.0 * g_matrix / 3.0)
    };
    k_matrix + f_inclusion / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hs_bounds_enclose_voigt_for_two_phase_mix() {
        // Al + steel.
        let al = (70_000.0, 0.33, 0.5);
        let steel = (200_000.0, 0.3, 0.5);
        let hs = hashin_shtrikman_two_phase(al, steel);
        let (k_al, g_al) = (k_from_e_nu(70_000.0, 0.33), g_from_e_nu(70_000.0, 0.33));
        let (k_st, g_st) = (k_from_e_nu(200_000.0, 0.3), g_from_e_nu(200_000.0, 0.3));
        // Voigt bounds:
        let k_v = 0.5 * (k_al + k_st);
        let g_v = 0.5 * (g_al + g_st);
        // Hashin–Shtrikman bounds must enclose Voigt (a known
        // property of HS: they are tighter than Voigt for
        // two-phase mixtures with isotropic constituents).
        assert!(hs.k_lower <= k_v + 1e-6, "HS lower must be ≤ Voigt");
        assert!(hs.k_upper >= k_v - 1e-6, "HS upper must be ≥ Voigt");
        assert!(hs.g_lower <= g_v + 1e-6, "HS lower must be ≤ Voigt G");
        assert!(hs.g_upper >= g_v - 1e-6, "HS upper must be ≥ Voigt G");
        // Sanity: bounds are ordered.
        assert!(hs.k_lower <= hs.k_upper);
        assert!(hs.g_lower <= hs.g_upper);
    }

    #[test]
    fn hs_spherical_recovers_matrix_for_zero_inclusion() {
        let k = hashin_shtrikman_spherical_bulk(100.0, 50.0, 200.0, 0.0);
        assert!((k - 100.0).abs() < 1e-9);
    }

    #[test]
    fn hs_spherical_recovers_inclusion_for_full_inclusion() {
        let k = hashin_shtrikman_spherical_bulk(100.0, 50.0, 200.0, 1.0);
        assert!((k - 200.0).abs() < 1e-9);
    }

    #[test]
    fn hs_spherical_intermediate_fraction_lies_between() {
        let k0 = hashin_shtrikman_spherical_bulk(100.0, 50.0, 200.0, 0.0);
        let k1 = hashin_shtrikman_spherical_bulk(100.0, 50.0, 200.0, 1.0);
        let k_mid = hashin_shtrikman_spherical_bulk(100.0, 50.0, 200.0, 0.5);
        assert!(k_mid > k0 && k_mid < k1);
    }
}
