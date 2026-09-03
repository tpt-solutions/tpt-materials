//! Rule-of-mixtures (engineering) estimate.
//!
//! A simple linear weighted average of the phase stiffnesses
//! `C_eff = Σ f_r C_r`.  Equivalent to the Voigt bound.

use tpt_mat_crystal_plasticity::SymmetricFourthOrder;

/// Rule-of-mixtures effective stiffness for an `N`-phase composite.
pub fn rule_of_mixtures(phases: &[(SymmetricFourthOrder, f64)]) -> SymmetricFourthOrder {
    let total_f: f64 = phases.iter().map(|(_, f)| *f).sum();
    if total_f <= 0.0 || phases.is_empty() {
        return SymmetricFourthOrder::new([[0.0_f64; 6]; 6]);
    }
    let mut c = [[0.0_f64; 6]; 6];
    for (c_r, f_r) in phases {
        let w = f_r / total_f;
        for i in 0..6 {
            for j in 0..6 {
                c[i][j] += w * c_r.data[i][j];
            }
        }
    }
    SymmetricFourthOrder::new(c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_phase_recovers_input() {
        let c = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let eff = rule_of_mixtures(&[(c.clone(), 1.0)]);
        for i in 0..6 {
            for j in 0..6 {
                assert!((eff.data[i][j] - c.data[i][j]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn two_phase_average_is_weighted() {
        let c_a = SymmetricFourthOrder::isotropic(70_000.0, 0.33);
        let c_b = SymmetricFourthOrder::isotropic(200_000.0, 0.3);
        let eff = rule_of_mixtures(&[(c_a.clone(), 0.3), (c_b.clone(), 0.7)]);
        let expected = 0.3 * c_a.data[0][0] + 0.7 * c_b.data[0][0];
        assert!((eff.data[0][0] - expected).abs() < 1e-6);
    }
}
