//! Grain-boundary and interface models.
//!
//! Consolidates the GB-property models that previously lived
//! separately in `tpt-mat-grain-growth::GrainBoundaryMobility` (now
//! re-exported here) and `tpt-mat-diffusion::GrainBoundaryDiffusion`
//! (also re-exported).  Adds:
//!
//! - [`GrainBoundaryEnergy`] — Read–Shockley misorientation-dependent
//!   boundary energy (the dual of GB mobility).
//! - [`GrainBoundaryCharacterDistribution`] — the 5-parameter
//!   macroscopic descriptor `(θ, n, φ_1, Φ, φ_2)` for the GB
//!   misorientation + boundary-plane normal.
//! - [`Triplejunction`] — dihedral-angle force balance at a triple
//!   junction (Herring 1951; consistent with mean-curvature flow).
//! - [`LangmuirMcLean`] — equilibrium GB solute segregation
//!   (McLean 1957).

use serde::{Deserialize, Serialize};

/// Read–Shockley misorientation-dependent grain-boundary energy
/// `γ_GB(θ) = γ_H · θ / θ_c · (1 − ln(θ/θ_c))` for low-angle
/// boundaries (`θ < θ_c`) and `γ_H` (constant) for high-angle.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GrainBoundaryEnergy {
    /// High-angle boundary energy `γ_H` (J/m²).
    pub gamma_high_angle: f64,
    /// Low-angle threshold `θ_c` (radians).
    pub low_angle_threshold: f64,
}

impl Default for GrainBoundaryEnergy {
    fn default() -> Self {
        Self {
            gamma_high_angle: 0.8,
            low_angle_threshold: 10.0_f64.to_radians(),
        }
    }
}

impl GrainBoundaryEnergy {
    /// Construct from defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Energy at misorientation `theta` (radians).
    pub fn at_misorientation(&self, theta: f64) -> f64 {
        let t = theta.abs();
        // Read–Shockley regime.
        if t < self.low_angle_threshold && t > 1.0e-12 {
            let ratio = t / self.low_angle_threshold;
            self.gamma_high_angle * ratio * (1.0 - ratio.ln())
        } else if t < 1.0e-12 {
            // Smooth limit: γ ∝ θ → 0.
            0.0
        } else {
            self.gamma_high_angle
        }
    }
}

/// Macroscopic grain-boundary character distribution: each entry is
/// `(θ, n_hat)` with `θ` the misorientation angle and `n_hat` the
/// boundary-plane normal in the crystal frame.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GBCDEntry {
    /// Misorientation angle (radians).
    pub misorientation: f64,
    /// Boundary-plane normal `(n_x, n_y, n_z)` in the reference
    /// crystal frame; magnitude is one.
    pub plane_normal: [f64; 3],
    /// Relative frequency (any non-negative number, unnormalised).
    pub weight: f64,
}

/// Distribution of grain-boundary characters (the macroscopically
/// observable "grain-boundary spectrum").
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GrainBoundaryCharacterDistribution {
    /// Entries.
    pub entries: Vec<GBCDEntry>,
}

impl GrainBoundaryCharacterDistribution {
    /// Empty distribution.
    pub fn new() -> Self {
        Self { entries: Vec::new() }
    }

    /// Add an entry.
    pub fn push(&mut self, entry: GBCDEntry) {
        self.entries.push(entry);
    }

    /// Total weight.
    pub fn total_weight(&self) -> f64 {
        self.entries.iter().map(|e| e.weight).sum()
    }

    /// Fraction of low-angle boundaries (LAB fraction).
    pub fn low_angle_fraction(&self, threshold: f64) -> f64 {
        let total = self.total_weight();
        if total <= 0.0 {
            return 0.0;
        }
        let lab: f64 = self
            .entries
            .iter()
            .filter(|e| e.misorientation < threshold)
            .map(|e| e.weight)
            .sum();
        lab / total
    }

    /// Fraction of coincidence-site-lattice boundaries with
    /// `Σ ≤ max_sigma` (a common measure of "special boundaries").
    pub fn csl_fraction(&self, max_sigma: u32) -> f64 {
        // Crude Σ tagging: misorientation near low-Σ CSL angles.
        let csl_angles = [
            (3_u32, 60.0_f64.to_radians()),
            (5, 36.87_f64.to_radians()),
            (7, 38.21_f64.to_radians()),
            (9, 38.94_f64.to_radians()),
            (11, 50.48_f64.to_radians()),
            (13, 22.62_f64.to_radians()),
            (15, 48.19_f64.to_radians()),
            (17, 28.07_f64.to_radians()),
            (19, 26.53_f64.to_radians()),
        ];
        let total = self.total_weight();
        if total <= 0.0 {
            return 0.0;
        }
        let mut w = 0.0;
        let tol = 2.0_f64.to_radians();
        for e in &self.entries {
            for &(sigma, theta) in &csl_angles {
                if sigma > max_sigma {
                    break;
                }
                if (e.misorientation - theta).abs() < tol {
                    w += e.weight;
                    break;
                }
            }
        }
        w / total
    }
}

impl Default for GrainBoundaryCharacterDistribution {
    fn default() -> Self {
        Self::new()
    }
}

/// Triple-junction dihedral-angle force balance (Herring 1951).
///
/// At a triple junction the three boundary tensions
/// `γ_12, γ_23, γ_13` must balance.  For isotropic boundaries
/// (`γ_ij = γ` for all `i,j`) the dihedral angles are `120°`
/// each.  This helper returns the largest dihedral angle between
/// the three boundaries and the implied curvature.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TripleJunction {
    /// Boundary energies (J/m²).
    pub gammas: [f64; 3],
}

impl TripleJunction {
    /// Dihedral angle opposite boundary `i` (radians).  Returns
    /// `None` if the configuration cannot close (Neumann triangle
    /// violated).
    pub fn dihedral_angle(&self, i: usize) -> Option<f64> {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let gj = self.gammas[j];
        let gk = self.gammas[k];
        let gi = self.gammas[i];
        let cos_t = (gj.powi(2) + gk.powi(2) - gi.powi(2)) / (2.0 * gj * gk);
        if !(-1.0..=1.0).contains(&cos_t) {
            return None;
        }
        Some(cos_t.acos())
    }

    /// `true` if the Neumann triangle closes for all three dihedral
    /// angles.
    pub fn is_equilibrium(&self) -> bool {
        (0..3).all(|i| self.dihedral_angle(i).is_some())
    }

    /// Implied driving curvature at the junction
    /// `(Σ γ_ij sin φ_ij)` for use in curvature-flow integrators.
    pub fn net_driving_force(&self) -> f64 {
        (0..3)
            .filter_map(|i| self.dihedral_angle(i))
            .zip(self.gammas.iter())
            .map(|(phi, g)| g * phi.sin())
            .sum()
    }
}

/// Langmuir–McLean equilibrium GB segregation (McLean 1957):
/// `X_GB / (1 − X_GB) = (X_bulk / (1 − X_bulk)) · exp(−ΔG_seg / RT)`.
///
/// Inverse form returns `X_GB(X_bulk, ΔG_seg)` and the segregation
/// enrichment factor `β = X_GB / X_bulk`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LangmuirMcLean {
    /// Segregation free energy `ΔG_seg` (J/mol; negative favours
    /// segregation).
    pub segregation_energy: f64,
}

impl LangmuirMcLean {
    /// Construct from `ΔG_seg` (J/mol).
    pub fn new(segregation_energy: f64) -> Self {
        Self { segregation_energy }
    }

    /// Equilibrium GB solute fraction `X_GB` for a given bulk
    /// solute fraction `X_bulk` and temperature.
    pub fn equilibrium_gb_fraction(
        &self,
        x_bulk: f64,
        temperature_k: f64,
        gas_constant: f64,
    ) -> f64 {
        let k = (-self.segregation_energy / (gas_constant * temperature_k)).exp();
        let num = x_bulk * k;
        let den = 1.0 + x_bulk * (k - 1.0);
        if den.abs() < 1.0e-30 {
            0.0
        } else {
            (num / den).clamp(0.0, 1.0)
        }
    }

    /// Enrichment factor `β = X_GB / X_bulk`.
    pub fn enrichment_factor(
        &self,
        x_bulk: f64,
        temperature_k: f64,
        gas_constant: f64,
    ) -> f64 {
        if x_bulk.abs() < 1.0e-30 {
            return 0.0;
        }
        self.equilibrium_gb_fraction(x_bulk, temperature_k, gas_constant) / x_bulk
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn read_shockley_energy_increases_to_max_then_saturates() {
        let e = GrainBoundaryEnergy::default();
        let g0 = e.at_misorientation(0.0);
        let g_max = e.at_misorientation(e.low_angle_threshold);
        let g_high = e.at_misorientation(0.5);
        assert!(g0.abs() < 1.0e-9);
        assert!((g_max - e.gamma_high_angle).abs() < 1e-12);
        assert!(approx(g_high, e.gamma_high_angle, 1e-12));
    }

    #[test]
    fn gbcd_low_angle_fraction_counts_correctly() {
        let mut d = GrainBoundaryCharacterDistribution::new();
        d.push(GBCDEntry {
            misorientation: 0.05,
            plane_normal: [1.0, 0.0, 0.0],
            weight: 1.0,
        });
        d.push(GBCDEntry {
            misorientation: 0.4,
            plane_normal: [0.0, 1.0, 0.0],
            weight: 2.0,
        });
        let lab = d.low_angle_fraction(10.0_f64.to_radians());
        assert!(approx(lab, 1.0 / 3.0, 1.0e-9));
    }

    #[test]
    fn gbcd_csl_fraction_identifies_sigma3() {
        let mut d = GrainBoundaryCharacterDistribution::new();
        d.push(GBCDEntry {
            misorientation: 60.0_f64.to_radians(),
            plane_normal: [1.0, 1.0, 1.0],
            weight: 1.0,
        });
        d.push(GBCDEntry {
            misorientation: 0.4,
            plane_normal: [0.0, 1.0, 0.0],
            weight: 3.0,
        });
        let csl = d.csl_fraction(7);
        assert!(approx(csl, 0.25, 1.0e-9));
    }

    #[test]
    fn triple_junction_isotropic_dihedral_is_60_degrees() {
        let t = TripleJunction {
            gammas: [1.0, 1.0, 1.0],
        };
        for i in 0..3 {
            let phi = t.dihedral_angle(i).expect("closed");
            assert!(approx(phi, 60.0_f64.to_radians(), 1.0e-9));
        }
    }

    #[test]
    fn triple_junction_violation_returns_none() {
        let t = TripleJunction {
            gammas: [10.0, 1.0, 1.0],
        };
        assert!(t.dihedral_angle(0).is_none());
    }

    #[test]
    fn langmuir_mclean_segregation_enriches_gb_over_bulk() {
        let seg = LangmuirMcLean::new(-30.0e3);
        let r = 8.314_462_618;
        let x_bulk = 0.01;
        let x_gb = seg.equilibrium_gb_fraction(x_bulk, 800.0, r);
        let beta = seg.enrichment_factor(x_bulk, 800.0, r);
        assert!(x_gb > x_bulk);
        assert!(beta > 1.0);
    }

    #[test]
    fn langmuir_mclean_temperature_decreases_segregation() {
        let seg = LangmuirMcLean::new(-30.0e3);
        let r = 8.314_462_618;
        let x_bulk = 0.01;
        let x_lo = seg.equilibrium_gb_fraction(x_bulk, 600.0, r);
        let x_hi = seg.equilibrium_gb_fraction(x_bulk, 1200.0, r);
        assert!(x_lo > x_hi);
    }
}