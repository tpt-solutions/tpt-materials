//! Effective diffusivity in porous / multiphase media.
//!
//! - [`effective_diffusivity`] — tortuosity-corrected
//!   `D_eff = φ D / τ²` where `τ` is the tortuosity.
//! - [`bruggeman_diffusivity`] — implicit Bruggeman relation
//!   `(1 − φ) (D_m − D_eff) / (D_m + (n − 1) D_eff) + φ (D_i − D_eff) / (D_i + (n − 1) D_eff) = 0`
//!   where `n` is the dimensionality (commonly 3 for an
//!   isotropic continuum).

/// Effective diffusivity in a porous medium of porosity `φ`
/// (volume fraction of voids) and intrinsic diffusivity `d_matrix`,
/// using a tortuosity `τ ≥ 1`:
///
/// `D_eff = φ D_matrix / τ²`
pub fn effective_diffusivity(porosity: f64, d_matrix: f64, tortuosity: f64) -> f64 {
    if porosity <= 0.0 {
        return 0.0;
    }
    let tau_sq = tortuosity * tortuosity;
    if tau_sq <= 0.0 {
        return 0.0;
    }
    porosity * d_matrix / tau_sq
}

/// Solve the Bruggeman implicit relation for `D_eff` via bisection.
///
/// `n` is the dimensionality (commonly 3).
pub fn bruggeman_diffusivity(porosity: f64, d_matrix: f64, d_inclusion: f64, n: f64) -> f64 {
    let f = porosity;
    let mut lo = d_matrix.min(d_inclusion);
    let mut hi = d_matrix.max(d_inclusion);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let lhs = (1.0 - f) * (d_matrix - mid) / (d_matrix + (n - 1.0) * mid)
            + f * (d_inclusion - mid) / (d_inclusion + (n - 1.0) * mid);
        if lhs > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn effective_diffusivity_zero_at_zero_porosity() {
        assert!(approx(
            effective_diffusivity(0.0, 1.0e-9, 1.5),
            0.0,
            1.0e-30
        ));
    }

    #[test]
    fn effective_diffusivity_decreases_with_tortuosity() {
        let d1 = effective_diffusivity(0.3, 1.0e-9, 1.0);
        let d2 = effective_diffusivity(0.3, 1.0e-9, 2.0);
        assert!(d2 < d1);
    }

    #[test]
    fn bruggeman_recovers_matrix_at_zero_porosity() {
        let d = bruggeman_diffusivity(0.0, 1.0e-9, 1.0e-12, 3.0);
        assert!(approx(d, 1.0e-9, 1.0e-15));
    }

    #[test]
    fn bruggeman_recovers_inclusion_at_full_porosity() {
        let d = bruggeman_diffusivity(1.0, 1.0e-9, 1.0e-12, 3.0);
        assert!(approx(d, 1.0e-12, 1.0e-15));
    }
}
