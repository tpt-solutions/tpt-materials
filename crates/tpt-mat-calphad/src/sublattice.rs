//! Sublattice model (Muggianu-style).
//!
//! `s` sublattices, each with `n_i` species of site fraction `y_{ij}`
//! such that `Σ_j y_{ij} = 1`.  Ideal mixing configurational entropy per
//! mole of atoms:
//!
//! `^id S^cfg = −R Σ_i a_i Σ_j y_{ij} ln y_{ij}`
//!
//! where `a_i` is the site multiplicity (atoms per formula unit on
//! sublattice `i`).  Configurational Gibbs energy contribution is
//! `−T · ^id S^cfg`.  End-member Gibbs energies per formula unit are
//! supplied by the caller.

use serde::{Deserialize, Serialize};

/// Single species on a sublattice.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SublatticeSpecies {
    /// Human-readable species name.
    pub name: String,
    /// Site fraction of this species on its sublattice.
    pub site_fraction: f64,
}

/// Sublattice definition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Sublattice {
    /// Site multiplicity (atoms per formula unit).
    pub multiplicity: f64,
    /// Species occupying this sublattice.
    pub species: Vec<SublatticeSpecies>,
}

/// Top-level sublattice configuration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SublatticeConfig {
    /// Sublattices in formula-unit order.
    pub sublattices: Vec<Sublattice>,
}

/// Muggianu-style sublattice model.
///
/// Holds the sublattice configuration and the temperature-dependent
/// end-member Gibbs energies of the *current* composition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SublatticeModel {
    /// Configuration.
    pub config: SublatticeConfig,
    /// Universal gas constant.
    pub r: f64,
    /// Temperature (K).
    pub t_k: f64,
}

impl Default for SublatticeModel {
    fn default() -> Self {
        // Example: FCC_A1 (Cu, Zn)_1 (Zn, Va)_1 — a hypothetical
        // brass-like two-sublattice model.
        Self {
            config: SublatticeConfig {
                sublattices: vec![
                    Sublattice {
                        multiplicity: 1.0,
                        species: vec![
                            SublatticeSpecies {
                                name: "Cu".into(),
                                site_fraction: 0.5,
                            },
                            SublatticeSpecies {
                                name: "Zn".into(),
                                site_fraction: 0.5,
                            },
                        ],
                    },
                    Sublattice {
                        multiplicity: 1.0,
                        species: vec![
                            SublatticeSpecies {
                                name: "Zn".into(),
                                site_fraction: 0.7,
                            },
                            SublatticeSpecies {
                                name: "Va".into(),
                                site_fraction: 0.3,
                            },
                        ],
                    },
                ],
            },
            r: 8.314_462_618,
            t_k: 1000.0,
        }
    }
}

impl SublatticeModel {
    /// Construct an explicit model.
    pub fn new(config: SublatticeConfig, t_k: f64) -> Self {
        Self {
            config,
            r: 8.314_462_618,
            t_k,
        }
    }

    /// Normalise all site fractions on each sublattice to sum to 1
    /// (in place).  Returns true if any normalisation was actually
    /// required.
    pub fn normalise(&mut self) -> bool {
        let mut changed = false;
        for sl in &mut self.config.sublattices {
            let sum: f64 = sl.species.iter().map(|s| s.site_fraction).sum();
            if (sum - 1.0).abs() > 1e-12 {
                changed = true;
                for s in &mut sl.species {
                    s.site_fraction /= sum.max(1e-12);
                }
            }
        }
        changed
    }

    /// Muggianu ideal-mixing configurational Gibbs energy per formula
    /// unit: `−T · ^id S^cfg`.
    pub fn configurational_gibbs(&self) -> f64 {
        let mut s_cfg = 0.0;
        for sl in &self.config.sublattices {
            for sp in &sl.species {
                let y = sp.site_fraction.clamp(1e-12, 1.0);
                s_cfg -= sl.multiplicity * y * y.ln();
            }
        }
        s_cfg *= self.r;
        -self.t_k * s_cfg
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalisation_sums_to_one() {
        let mut m = SublatticeModel::default();
        // Perturb a site fraction.
        m.config.sublattices[0].species[0].site_fraction = 0.7;
        m.config.sublattices[0].species[1].site_fraction = 0.4;
        m.normalise();
        for sl in &m.config.sublattices {
            let s: f64 = sl.species.iter().map(|s| s.site_fraction).sum();
            assert!((s - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn configurational_gibbs_is_nonpositive() {
        // Ideal mixing always lowers G (or leaves it unchanged at the
        // endpoints); for mixed fractions it must be strictly negative.
        let m = SublatticeModel::default();
        assert!(m.configurational_gibbs() <= 0.0);
    }

    #[test]
    fn pure_endpoints_have_zero_gibbs() {
        let config = SublatticeConfig {
            sublattices: vec![Sublattice {
                multiplicity: 1.0,
                species: vec![SublatticeSpecies {
                    name: "Cu".into(),
                    site_fraction: 1.0,
                }],
            }],
        };
        let m = SublatticeModel::new(config, 1000.0);
        assert_eq!(m.configurational_gibbs(), 0.0);
    }
}
