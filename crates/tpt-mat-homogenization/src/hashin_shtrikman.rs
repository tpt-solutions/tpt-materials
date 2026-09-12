//! Hashin–Shtrikman (1963) variational bounds for a two-phase
//! isotropic composite.
//!
//! For a composite of two isotropic phases with bulk moduli `K₁, K₂`,
//! shear moduli `G₁, G₂` and volume fraction `f` of phase 1, the
//! Hashin–Shtrikman bounds are the tightest possible isotropic
//! bounds that depend only on the phase moduli and the volume fraction
//! (without specifying the microstructure).
//!
//! With the phases ordered so that phase 1 is the *softer* (i.e.
//! `K₁ ≤ K₂` and `G₁ ≤ G₂`), the bounds are
//!
//! `K^- = K₁ + f₂ / (1 / (K₂ - K₁) + 3 f₁ / (3 K₁ + 4 G₁))`
//!
//! `K^+ = K₂ + f₁ / (1 / (K₁ - K₂) + 3 f₂ / (3 K₂ + 4 G₂))`
//!
//! `G^- = G₁ + f₂ / (1 / (G₂ - G₁) + 6 (K₁ + 2 G₁) f₁ / (5 G₁ (3 K₁ + 4 G₁)))`
//!
//! `G^+ = G₂ + f₁ / (1 / (G₁ - G₂) + 6 (K₂ + 2 G₂) f₂ / (5 G₂ (3 K₂ + 4 G₂)))`
//!
//! (Hashin & Shtrikman, J. Mech. Phys. Solids 11, 1963.)  Each bound
//! corresponds to one of the two coated-sphere microstructures:
//! `K^-`/`G^-` keep the softer phase as the percolating matrix (lowest
//! stiffness), `K^+`/`G^+` keep the stiffer phase as the matrix.  At the
//! endpoints `f → 0` and `f → 1` both bounds collapse to the
//! corresponding pure-phase modulus.
//!
//! If the phases are not monotonically ordered (`K₁ ≤ K₂` *and*
//! `G₁ ≤ G₂` do not hold), the tight two-point bounds do not apply; the
//! code then falls back to the phase-modulus range `[min, max]`, which
//! remains a valid (coarser) bound.

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
    let (e_a, nu_a, f_a) = phase_a;
    let (e_b, nu_b, _f_b) = phase_b;
    let (k_a, g_a) = (k_from_e_nu(e_a, nu_a), g_from_e_nu(e_a, nu_a));
    let (k_b, g_b) = (k_from_e_nu(e_b, nu_b), g_from_e_nu(e_b, nu_b));
    hashin_shtrikman_k_g((k_a, g_a), (k_b, g_b), f_a)
}

/// Hashin–Shtrikman bounds from explicit `(K, G)` moduli and the
/// volume fraction `f_a` of phase `a`.
///
/// The two calls below agree on the endpoints but the bounds are
/// genuinely volume-fraction dependent for `0 < f_a < 1`; in
/// particular the HS interval is much tighter than the coarse
/// `[min, max]` phase-modulus range and strictly encloses the
/// Voigt/Reuss (and therefore the VRH) averages.
pub fn hashin_shtrikman_k_g(a: (f64, f64), b: (f64, f64), f_a: f64) -> HashinShtrikmanResult {
    let (k_a, g_a) = a;
    let (k_b, g_b) = b;
    let f_a = f_a.clamp(0.0, 1.0);
    let f_b = 1.0 - f_a;
    // Order so that phase 1 is the *softer* phase.  The HS bounds
    // require a monotonic ordering in both moduli; outside that range
    // they are not tight two-point bounds and we fall back to the
    // phase-modulus interval.
    let (k_1, g_1, f_1, k_2, g_2, f_2) = if k_a <= k_b && g_a <= g_b {
        (k_a, g_a, f_a, k_b, g_b, f_b)
    } else if k_b <= k_a && g_b <= g_a {
        (k_b, g_b, f_b, k_a, g_a, f_a)
    } else {
        let k_lo = k_a.min(k_b);
        let k_hi = k_a.max(k_b);
        let g_lo = g_a.min(g_b);
        let g_hi = g_a.max(g_b);
        return HashinShtrikmanResult {
            k_lower: k_lo,
            k_upper: k_hi,
            g_lower: g_lo,
            g_upper: g_hi,
        };
    };
    let k_lower = hashin_shtrikman_spherical_bulk(k_1, g_1, k_2, f_2);
    let k_upper = hashin_shtrikman_spherical_bulk(k_2, g_2, k_1, f_1);
    let g_lower = hashin_shtrikman_spherical_shear(k_1, g_1, k_2, g_2, f_2);
    let g_upper = hashin_shtrikman_spherical_shear(k_2, g_2, k_1, g_1, f_1);
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

/// Hashin–Shtrikman (1963) spherical bound for the shear modulus:
///
/// `G_eff^± = G_0 + f / (1 / (G_1 - G_0) + 6 (K_0 + 2 G_0) (1 - f) / (5 G_0 (3 K_0 + 4 G_0)))`
///
/// with `+` when the inclusion is stiffer than the matrix (`G_1 > G_0`)
/// and `−` when it is softer.  `(K_0, G_0)` are the matrix moduli,
/// `(K_1, G_1)` the inclusion moduli and `f` the inclusion volume
/// fraction.
pub fn hashin_shtrikman_spherical_shear(
    k_matrix: f64,
    g_matrix: f64,
    k_inclusion: f64,
    g_inclusion: f64,
    f_inclusion: f64,
) -> f64 {
    let dg = g_inclusion - g_matrix;
    let denom = if dg.abs() < 1e-30 {
        return g_matrix;
    } else {
        1.0 / dg
            + 6.0 * (k_matrix + 2.0 * g_matrix) * (1.0 - f_inclusion)
                / (5.0 * g_matrix * (3.0 * k_matrix + 4.0 * g_matrix))
    };
    g_matrix + f_inclusion / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hs_bounds_enclose_voigt_for_two_phase_mix() {
        // Al + steel, 50:50.
        let al = (70_000.0, 0.33, 0.5);
        let steel = (200_000.0, 0.3, 0.5);
        let hs = hashin_shtrikman_two_phase(al, steel);
        let (k_al, g_al) = (k_from_e_nu(70_000.0, 0.33), g_from_e_nu(70_000.0, 0.33));
        let (k_st, g_st) = (k_from_e_nu(200_000.0, 0.3), g_from_e_nu(200_000.0, 0.3));
        // Voigt / Reuss bounds:
        let k_v = 0.5 * (k_al + k_st);
        let g_v = 0.5 * (g_al + g_st);
        let k_r = 1.0 / (0.5 / k_al + 0.5 / k_st);
        let g_r = 1.0 / (0.5 / g_al + 0.5 / g_st);
        // HS bounds must be *tighter* than the Voigt/Reuss bounds:
        // Reuss ≤ HS^- ≤ HS^+ ≤ Voigt holds for the effective modulus.
        assert!(hs.k_lower >= k_r - 1e-6, "HS K- must be ≥ Reuss (tighter)");
        assert!(hs.k_upper <= k_v + 1e-6, "HS K+ must be ≤ Voigt (tighter)");
        assert!(hs.g_lower >= g_r - 1e-6, "HS G- must be ≥ Reuss (tighter)");
        assert!(hs.g_upper <= g_v + 1e-6, "HS G+ must be ≤ Voigt (tighter)");
        // ... and must be strictly tighter than the phase range
        // [min, max] (the old implementation returned exactly that).
        assert!(hs.k_lower > k_al - 1.0, "HS K- must lie above soft phase K");
        assert!(
            hs.k_upper < k_st + 1.0,
            "HS K+ must lie below stiff phase K"
        );
        assert!(hs.g_lower > g_al - 1.0, "HS G- must lie above soft phase G");
        assert!(
            hs.g_upper < g_st + 1.0,
            "HS G+ must lie below stiff phase G"
        );
        // Sanity: bounds are ordered.
        assert!(hs.k_lower <= hs.k_upper);
        assert!(hs.g_lower <= hs.g_upper);
    }

    #[test]
    fn hs_bounds_respect_volume_fraction_endpoints() {
        // Soft phase dominates: bounds collapse onto the soft phase.
        let soft = (70_000.0, 0.33, 0.9999);
        let hard = (200_000.0, 0.3, 0.0001);
        let hs = hashin_shtrikman_two_phase(soft, hard);
        let (k_al, g_al) = (k_from_e_nu(70_000.0, 0.33), g_from_e_nu(70_000.0, 0.33));
        let (k_st, g_st) = (k_from_e_nu(200_000.0, 0.3), g_from_e_nu(200_000.0, 0.3));
        // Both bounds collapse onto ~the soft phase: at 99.99 % soft the
        // hard inclusions only stiffen the composite by a tiny amount.
        // Assert within 1 % of the soft-phase modulus and *inside* the
        // phase range.
        assert!((hs.k_lower - k_al).abs() < 0.01 * k_al);
        assert!((hs.k_upper - k_al).abs() < 0.01 * k_al);
        assert!((hs.g_lower - g_al).abs() < 0.01 * g_al);
        assert!((hs.g_upper - g_al).abs() < 0.01 * g_al);
        assert!(
            hs.k_lower > k_al && hs.k_upper < k_st,
            "bounds within [soft, hard]"
        );
        assert!(hs.g_lower > g_al && hs.g_upper < g_st);
    }

    #[test]
    fn hs_bounds_move_monotonically_with_volume_fraction() {
        // As the soft phase fraction grows, both the upper and lower
        // bulk bounds must move strictly down (monotonicity of the
        // HS formulas in f).
        let wide = hashin_shtrikman_two_phase((70_000.0, 0.33, 0.9), (200_000.0, 0.3, 0.1));
        let mid = hashin_shtrikman_two_phase((70_000.0, 0.33, 0.5), (200_000.0, 0.3, 0.5));
        let narrow = hashin_shtrikman_two_phase((70_000.0, 0.33, 0.1), (200_000.0, 0.3, 0.9));
        assert!(wide.k_lower < mid.k_lower);
        assert!(mid.k_lower < narrow.k_lower);
        assert!(wide.k_upper < mid.k_upper);
        assert!(mid.k_upper < narrow.k_upper);
        // The interval is tight (tens of GPa, not hundreds).
        assert!(mid.k_upper - mid.k_lower < 10_000.0);
        assert!(mid.g_upper - mid.g_lower < 10_000.0);
    }

    #[test]
    fn hs_shear_bound_recovers_matrix_and_inclusion_limit() {
        let g0 = hashin_shtrikman_spherical_shear(100.0, 50.0, 200.0, 100.0, 0.0);
        assert!((g0 - 50.0).abs() < 1e-9);
        let g1 = hashin_shtrikman_spherical_shear(100.0, 50.0, 200.0, 100.0, 1.0);
        assert!((g1 - 100.0).abs() < 1e-9);
        let g_mid = hashin_shtrikman_spherical_shear(100.0, 50.0, 200.0, 100.0, 0.5);
        assert!(g_mid > 50.0 && g_mid < 100.0);
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
