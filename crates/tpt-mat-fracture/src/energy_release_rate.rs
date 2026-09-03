//! Energy release rate (Irwin relation + J-integral helper).

use serde::{Deserialize, Serialize};

/// Plane-stress / plane-strain elastic modulus used in the Irwin
/// relation `G = K² / E'`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum IrwinModulus {
    /// Plane-stress modulus `E' = E`.
    PlaneStress {
        /// Young's modulus (Pa).
        youngs_modulus: f64,
    },
    /// Plane-strain modulus `E' = E / (1 − ν²)`.
    PlaneStrain {
        /// Young's modulus (Pa).
        youngs_modulus: f64,
        /// Poisson's ratio.
        poissons_ratio: f64,
    },
}

impl IrwinModulus {
    /// Reduced modulus `E'` (Pa).
    pub fn value(&self) -> f64 {
        match self {
            IrwinModulus::PlaneStress { youngs_modulus } => *youngs_modulus,
            IrwinModulus::PlaneStrain {
                youngs_modulus,
                poissons_ratio,
            } => youngs_modulus / (1.0 - poissons_ratio * poissons_ratio),
        }
    }
}

/// Irwin energy release rate `G = K_I² / E'`.
///
/// # Arguments
///
/// - `k_i`     — mode-I stress-intensity factor (Pa·m^0.5)
/// - `modulus` — reduced elastic modulus
pub fn energy_release_rate_irwin(k_i: f64, modulus: &IrwinModulus) -> f64 {
    k_i * k_i / modulus.value()
}

/// J-integral-based energy release rate proxy from a far-field
/// `J`-value (commonly obtained from a contour integral around
/// the crack tip).
///
/// This function simply returns the input, after a positivity
/// guard.  It exists so the cohesive-zone / phase-field code
/// paths can dispatch on the same `EnergyReleaseRate` enum.
///
/// # Arguments
///
/// - `j_value` — value of the J-integral computed around the
///   crack tip (J/m²)
pub fn energy_release_rate_j_integral(j_value: f64) -> f64 {
    j_value.max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn plane_stress_modulus_is_e() {
        let m = IrwinModulus::PlaneStress {
            youngs_modulus: 200.0e9,
        };
        assert!(approx(m.value(), 200.0e9, 1.0e-6));
    }

    #[test]
    fn plane_strain_modulus_is_e_over_one_minus_nu_sq() {
        let m = IrwinModulus::PlaneStrain {
            youngs_modulus: 200.0e9,
            poissons_ratio: 0.3,
        };
        let expected = 200.0e9 / (1.0 - 0.3 * 0.3);
        assert!(approx(m.value(), expected, 1.0e-3));
    }

    #[test]
    fn irwin_g_increases_with_k_squared() {
        let m = IrwinModulus::PlaneStress {
            youngs_modulus: 200.0e9,
        };
        let g1 = energy_release_rate_irwin(1.0e6, &m);
        let g2 = energy_release_rate_irwin(2.0e6, &m);
        assert!(approx(g2, 4.0 * g1, 1.0e-3));
    }

    #[test]
    fn irwin_g_lower_in_plane_strain() {
        // Larger E' -> smaller G.
        let plane_stress = IrwinModulus::PlaneStress {
            youngs_modulus: 200.0e9,
        };
        let plane_strain = IrwinModulus::PlaneStrain {
            youngs_modulus: 200.0e9,
            poissons_ratio: 0.3,
        };
        let g_ps = energy_release_rate_irwin(1.0e6, &plane_stress);
        let g_pe = energy_release_rate_irwin(1.0e6, &plane_strain);
        assert!(g_pe < g_ps);
    }
}