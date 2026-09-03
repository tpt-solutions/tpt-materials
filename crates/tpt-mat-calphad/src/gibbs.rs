//! Gibbs energy of a binary A–B phase: end-member + ideal-mixing +
//! Redlich–Kister excess.

use serde::{Deserialize, Serialize};

use crate::redlich_kister::{RedlichKister, RedlichKisterParams};

/// End-member Gibbs energy `G^A(T)` and `G^B(T)` at the temperature
/// of interest.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EndMember {
    /// `G^A(T)` (energy / mole).
    pub g_a: f64,
    /// `G^B(T)` (energy / mole).
    pub g_b: f64,
}

/// Ideal-mixing contribution:
///
/// `^id G(x_B) = R T (x_B ln x_B + (1 − x_B) ln (1 − x_B))`
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct IdealMixing {
    /// Universal gas constant `R` (energy / (mol·K)).
    pub r: f64,
    /// Temperature (K).
    pub t_k: f64,
}

impl Default for IdealMixing {
    fn default() -> Self {
        Self {
            r: 8.314_462_618,
            t_k: 1000.0,
        }
    }
}

impl IdealMixing {
    /// `^id G(x_B)` (returns 0 for `x_B ∈ {0, 1}`).
    pub fn contribution(&self, x_b: f64) -> f64 {
        if x_b <= 0.0 || x_b >= 1.0 {
            return 0.0;
        }
        self.r * self.t_k * (x_b * x_b.ln() + (1.0 - x_b) * (1.0 - x_b).ln())
    }

    /// `d^id G / dx_B = R T (ln x_B − ln (1 − x_B))`.
    pub fn derivative(&self, x_b: f64) -> f64 {
        if x_b <= 0.0 || x_b >= 1.0 {
            return 0.0;
        }
        self.r * self.t_k * (x_b.ln() - (1.0 - x_b).ln())
    }
}

/// Gibbs-energy model for a binary phase.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GibbsEnergyModel {
    /// End-member Gibbs energies at the temperature of interest.
    pub end_member: EndMember,
    /// Ideal-mixing contribution.
    pub ideal: IdealMixing,
    /// Redlich–Kister excess.
    pub excess: RedlichKister,
}

impl Default for GibbsEnergyModel {
    fn default() -> Self {
        Self {
            end_member: EndMember {
                g_a: -30000.0,
                g_b: -25000.0,
            },
            ideal: IdealMixing::default(),
            excess: RedlichKister::default(),
        }
    }
}

impl GibbsEnergyModel {
    /// End-member linear interpolation `x_A G^A + x_B G^B`.
    pub fn reference(&self, x_b: f64) -> f64 {
        let xa = 1.0 - x_b;
        xa * self.end_member.g_a + x_b * self.end_member.g_b
    }

    /// Total Gibbs energy `G(x_B) = G_ref + ^id G + ^E G`.
    pub fn total(&self, x_b: f64) -> f64 {
        self.reference(x_b) + self.ideal.contribution(x_b) + self.excess.excess(x_b)
    }

    /// `dG / dx_B`.
    pub fn derivative(&self, x_b: f64) -> f64 {
        let d_ref = self.end_member.g_b - self.end_member.g_a;
        d_ref + self.ideal.derivative(x_b) + self.excess.excess_derivative(x_b)
    }

    /// Re-build the model at a new temperature (replacing the
    /// end-member values and the ideal-mixing temperature).
    pub fn at_temperature(end_member: EndMember, t_k: f64, excess: RedlichKisterParams) -> Self {
        Self {
            end_member,
            ideal: IdealMixing {
                r: 8.314_462_618,
                t_k,
            },
            excess: RedlichKister::new(excess),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ideal_mixing_is_zero_at_endpoints_and_negative_inside() {
        let im = IdealMixing::default();
        assert_eq!(im.contribution(0.0), 0.0);
        assert_eq!(im.contribution(1.0), 0.0);
        let mid = im.contribution(0.5);
        assert!(mid < 0.0);
    }

    #[test]
    fn total_gibbs_is_reference_plus_extras() {
        let m = GibbsEnergyModel::default();
        for &x in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let total = m.total(x);
            let ref_ = m.reference(x);
            let id = m.ideal.contribution(x);
            let ex = m.excess.excess(x);
            assert!((total - (ref_ + id + ex)).abs() < 1e-12);
        }
    }
}
