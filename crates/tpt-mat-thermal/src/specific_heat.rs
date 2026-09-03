//! Effective specific heat: mass-weighted rule of mixtures.

/// Mass-weighted rule-of-mixtures specific heat.
///
/// `cp_eff = Σ (f_i ρ_i cp_i) / Σ (f_i ρ_i)`
/// where `f_i` is volume fraction, `ρ_i` density, and
/// `cp_i` specific heat.
///
/// For `n = 1`, simply returns `cp_phases[0]`.
pub fn effective_specific_heat(cp_phases: &[f64], densities: &[f64], f: &[f64]) -> f64 {
    assert_eq!(cp_phases.len(), f.len(), "cp_phases and f must match");
    assert_eq!(densities.len(), f.len(), "densities and f must match");
    let num: f64 = f
        .iter()
        .zip(cp_phases.iter().zip(densities.iter()))
        .map(|(fi, (cp_i, rho_i))| fi * cp_i * rho_i)
        .sum();
    let denom: f64 = f.iter().zip(densities.iter()).map(|(fi, rho_i)| fi * rho_i).sum();
    if denom.abs() < 1.0e-30 {
        return 0.0;
    }
    num / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn cp_reduces_to_single_phase_value() {
        let cp = effective_specific_heat(&[500.0], &[8000.0], &[1.0]);
        assert!(approx(cp, 500.0, 1.0e-9));
    }

    #[test]
    fn cp_mass_weighted_two_phases() {
        // 50% volume A: cp=500, ρ=8000; 50% B: cp=900, ρ=2700.
        // Mass fractions: m_A/m = 0.5*8000/(0.5*8000+0.5*2700) = 0.748
        // cp_eff = 0.748*500 + 0.252*900 = 374 + 227 = 601
        let cp = effective_specific_heat(&[500.0, 900.0], &[8000.0, 2700.0], &[0.5, 0.5]);
        let m_a = 0.5 * 8000.0;
        let m_b = 0.5 * 2700.0;
        let expected = (m_a * 500.0 + m_b * 900.0) / (m_a + m_b);
        assert!(approx(cp, expected, 1.0e-9));
    }
}