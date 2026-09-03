//! Bulk free-energy densities and total free-energy functionals.

use serde::{Deserialize, Serialize};

/// Bulk free-energy density `f(η)` for Allen-Cahn models.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum BulkEnergy {
    /// Symmetric double-well `f(η) = W η² (1 − η)²`.
    DoubleWell {
        /// Well depth `W`.
        well_depth: f64,
    },
    /// 6th-order polynomial `f(η) = W (η² (1 − η)² + a (η − 1/2))`.
    Polynomial {
        /// Well depth `W`.
        well_depth: f64,
        /// Asymmetry coefficient.
        asymmetry: f64,
    },
    /// Regular-solution `f(c) = Ω c (1 − c) + RT [c ln c + (1 − c) ln(1 − c)]`.
    RegularSolution(RegularSolutionParams),
}

impl BulkEnergy {
    /// Evaluate `f(x)`.
    pub fn value(&self, x: f64) -> f64 {
        match *self {
            BulkEnergy::DoubleWell { well_depth: w } => w * x * x * (1.0 - x).powi(2),
            BulkEnergy::Polynomial {
                well_depth: w,
                asymmetry: a,
            } => w * (x * x * (1.0 - x).powi(2) + a * (x - 0.5)),
            BulkEnergy::RegularSolution(p) => p.value(x),
        }
    }

    /// First derivative `f′(x)`.
    pub fn derivative(&self, x: f64) -> f64 {
        match *self {
            BulkEnergy::DoubleWell { well_depth: w } => {
                let u = x * (1.0 - x);
                2.0 * w * u * (1.0 - 2.0 * x)
            }
            BulkEnergy::Polynomial {
                well_depth: w,
                asymmetry: a,
            } => {
                let u = x * (1.0 - x);
                2.0 * w * u * (1.0 - 2.0 * x) + w * a
            }
            BulkEnergy::RegularSolution(p) => p.derivative(x),
        }
    }

    /// Second derivative `f″(x)`.
    pub fn second_derivative(&self, x: f64) -> f64 {
        match *self {
            BulkEnergy::DoubleWell { well_depth: w } => {
                let u = x * (1.0 - x);
                2.0 * w * u * (-1.0 - 4.0 * x)
            }
            BulkEnergy::Polynomial {
                well_depth: w,
                asymmetry: _,
            } => {
                let u = x * (1.0 - x);
                2.0 * w * u * (-1.0 - 4.0 * x)
            }
            BulkEnergy::RegularSolution(p) => p.second_derivative(x),
        }
    }
}

impl Default for BulkEnergy {
    fn default() -> Self {
        BulkEnergy::DoubleWell { well_depth: 1.0 }
    }
}

/// Regular-solution free-energy parameters.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RegularSolutionParams {
    /// Interaction parameter `Ω`.
    pub omega: f64,
    /// Gas constant `RT` (energy / mole).
    pub rt: f64,
}

impl RegularSolutionParams {
    /// Evaluate `f(c) = Ω c (1 − c) + RT [c ln c + (1 − c) ln(1 − c)]`.
    pub fn value(&self, c: f64) -> f64 {
        let eps = 1.0e-12;
        let c = c.clamp(eps, 1.0 - eps);
        self.omega * c * (1.0 - c) + self.rt * (c * c.ln() + (1.0 - c) * (1.0 - c).ln())
    }

    /// `f′(c) = Ω (1 − 2c) + RT ln(c / (1 − c))`.
    pub fn derivative(&self, c: f64) -> f64 {
        let eps = 1.0e-12;
        let c = c.clamp(eps, 1.0 - eps);
        self.omega * (1.0 - 2.0 * c) + self.rt * (c / (1.0 - c)).ln()
    }

    /// `f″(c) = −2 Ω + RT / (c (1 − c))`.
    pub fn second_derivative(&self, c: f64) -> f64 {
        let eps = 1.0e-12;
        let c = c.clamp(eps, 1.0 - eps);
        -2.0 * self.omega + self.rt / (c * (1.0 - c))
    }
}

/// Free-energy functional:
///
/// `F[η, c] = ∫ [f_bulk(η) + (κ/2) |∇η|² + f_chem(c)] dV`
///
/// evaluated on a regular 2D grid via central differences.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FreeEnergyFunctional {
    /// Gradient-energy coefficient `κ`.
    pub gradient_coefficient: f64,
    /// Cell spacing.
    pub dx: f64,
}

impl FreeEnergyFunctional {
    /// Compute the total free energy for a 2D field `eta` (length
    /// `nx * ny`, row-major) under the given bulk energy.
    pub fn evaluate_2d(&self, eta: &[f64], nx: usize, ny: usize, bulk: &BulkEnergy) -> f64 {
        let dx = self.dx;
        let mut total = 0.0;
        for i in 0..ny {
            for j in 0..nx {
                let c = i * nx + j;
                let e = eta[c];
                total += bulk.value(e);
                // Gradient energy: (κ/2) ( (∂x e)² + (∂y e)² ).
                let ex = if j + 1 < nx { eta[c + 1] } else { eta[c - 1] };
                let emx = if j > 0 { eta[c - 1] } else { eta[c + 1] };
                let ey = if i + 1 < ny { eta[c + nx] } else { eta[c - nx] };
                let emy = if i > 0 { eta[c - nx] } else { eta[c + nx] };
                let gx = (ex - emx) * 0.5 / dx;
                let gy = (ey - emy) * 0.5 / dx;
                total += 0.5 * self.gradient_coefficient * (gx * gx + gy * gy);
            }
        }
        total * dx * dx
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn double_well_minima_at_zero_and_one() {
        let b = BulkEnergy::DoubleWell { well_depth: 1.0 };
        assert!(b.value(0.0).abs() < 1e-12);
        assert!(b.value(1.0).abs() < 1e-12);
        assert!(b.value(0.5) > 0.0);
    }

    #[test]
    fn double_well_derivative_zero_at_minima() {
        let b = BulkEnergy::DoubleWell { well_depth: 1.0 };
        assert!(b.derivative(0.0).abs() < 1e-12);
        assert!(b.derivative(1.0).abs() < 1e-12);
        // f' has zeros at 0, 0.5, 1; check the non-zero behaviour:
        // f'(x) = 2W x(1-x)(1-2x), so for x in (0.5, 1), f'(x) < 0
        // (drives toward x = 0), and for x in (0, 0.5), f'(x) > 0
        // (drives toward x = 1).
        assert!(b.derivative(0.4) > 0.0);
        assert!(b.derivative(0.6) < 0.0);
    }

    #[test]
    fn regular_solution_is_symmetric_about_half() {
        let p = RegularSolutionParams {
            omega: 10.0,
            rt: 1.0,
        };
        let f_low = p.value(0.2);
        let f_high = p.value(0.8);
        assert!((f_low - f_high).abs() < 1e-9);
        assert!((p.derivative(0.2) + p.derivative(0.8)).abs() < 1e-9);
    }
}
