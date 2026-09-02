//! Minimal periodic table for `tpt-materials`.
//!
//! Provides `atomic_mass(symbol)` and `atomic_radius(symbol)` for the
//! elements commonly encountered in alloys (light + transition + rare
//! earth).  Values follow IUPAC 2021 standard atomic weights (rounded
//! to 4 significant figures) and metallic radii from
//! Cordero et al., Dalton Trans. 2008 (rounded to 2 significant
//! figures, in pm).

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// Per-element atomic data: standard atomic mass (kg/mol) and metallic
/// radius (metres).  `None` for radius indicates a noble gas or
/// non-metallic element where no consensus metallic radius exists.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AtomicData {
    /// Standard atomic mass in kg/mol.
    pub mass_kg_per_mol: f64,
    /// Metallic radius in metres (`None` for non-metals).
    pub metallic_radius_m: Option<f64>,
}

/// Static lookup for atomic data.
///
/// Phase 1 ships H, Li, Be, B, C, N, O, Mg, Al, Si, Ti, V, Cr, Mn, Fe,
/// Co, Ni, Cu, Zn, Y, Zr, Nb, Mo, Ag, Sn, Hf, Ta, W, Pb.  Phase 2+ will
/// expand to the full periodic table.
#[derive(Debug, Clone, Copy)]
pub struct PeriodicTable;

impl PeriodicTable {
    /// Atomic data for a given symbol (case-sensitive, e.g. `"Fe"`).
    pub fn data(symbol: &str) -> Option<AtomicData> {
        Self::table().get(symbol).copied()
    }

    /// Standard atomic mass in kg/mol.
    pub fn atomic_mass(symbol: &str) -> Option<f64> {
        Self::data(symbol).map(|d| d.mass_kg_per_mol)
    }

    /// Metallic radius in metres.
    pub fn atomic_radius(symbol: &str) -> Option<f64> {
        Self::data(symbol).and_then(|d| d.metallic_radius_m)
    }

    fn table() -> &'static HashMap<&'static str, AtomicData> {
        use AtomicData as A;
        static TABLE: std::sync::OnceLock<HashMap<&'static str, AtomicData>> =
            std::sync::OnceLock::new();
        TABLE.get_or_init(|| {
            let mut t = HashMap::new();
            // mass kg/mol (IUPAC 2021 rounded), radius m (Cordero 2008 rounded)
            let rows: &[(&str, A)] = &[
                (
                    "H",
                    A {
                        mass_kg_per_mol: 1.008e-3,
                        metallic_radius_m: Some(25e-12),
                    },
                ),
                (
                    "Li",
                    A {
                        mass_kg_per_mol: 6.94e-3,
                        metallic_radius_m: Some(152e-12),
                    },
                ),
                (
                    "Be",
                    A {
                        mass_kg_per_mol: 9.0122e-3,
                        metallic_radius_m: Some(112e-12),
                    },
                ),
                (
                    "B",
                    A {
                        mass_kg_per_mol: 10.81e-3,
                        metallic_radius_m: None,
                    },
                ),
                (
                    "C",
                    A {
                        mass_kg_per_mol: 12.011e-3,
                        metallic_radius_m: None,
                    },
                ),
                (
                    "N",
                    A {
                        mass_kg_per_mol: 14.007e-3,
                        metallic_radius_m: None,
                    },
                ),
                (
                    "O",
                    A {
                        mass_kg_per_mol: 15.999e-3,
                        metallic_radius_m: None,
                    },
                ),
                (
                    "Mg",
                    A {
                        mass_kg_per_mol: 24.305e-3,
                        metallic_radius_m: Some(160e-12),
                    },
                ),
                (
                    "Al",
                    A {
                        mass_kg_per_mol: 26.982e-3,
                        metallic_radius_m: Some(143e-12),
                    },
                ),
                (
                    "Si",
                    A {
                        mass_kg_per_mol: 28.085e-3,
                        metallic_radius_m: None,
                    },
                ),
                (
                    "Ti",
                    A {
                        mass_kg_per_mol: 47.867e-3,
                        metallic_radius_m: Some(146e-12),
                    },
                ),
                (
                    "V",
                    A {
                        mass_kg_per_mol: 50.942e-3,
                        metallic_radius_m: Some(132e-12),
                    },
                ),
                (
                    "Cr",
                    A {
                        mass_kg_per_mol: 51.996e-3,
                        metallic_radius_m: Some(128e-12),
                    },
                ),
                (
                    "Mn",
                    A {
                        mass_kg_per_mol: 54.938e-3,
                        metallic_radius_m: Some(127e-12),
                    },
                ),
                (
                    "Fe",
                    A {
                        mass_kg_per_mol: 55.845e-3,
                        metallic_radius_m: Some(126e-12),
                    },
                ),
                (
                    "Co",
                    A {
                        mass_kg_per_mol: 58.933e-3,
                        metallic_radius_m: Some(125e-12),
                    },
                ),
                (
                    "Ni",
                    A {
                        mass_kg_per_mol: 58.693e-3,
                        metallic_radius_m: Some(124e-12),
                    },
                ),
                (
                    "Cu",
                    A {
                        mass_kg_per_mol: 63.546e-3,
                        metallic_radius_m: Some(128e-12),
                    },
                ),
                (
                    "Zn",
                    A {
                        mass_kg_per_mol: 65.38e-3,
                        metallic_radius_m: Some(134e-12),
                    },
                ),
                (
                    "Y",
                    A {
                        mass_kg_per_mol: 88.906e-3,
                        metallic_radius_m: Some(180e-12),
                    },
                ),
                (
                    "Zr",
                    A {
                        mass_kg_per_mol: 91.224e-3,
                        metallic_radius_m: Some(160e-12),
                    },
                ),
                (
                    "Nb",
                    A {
                        mass_kg_per_mol: 92.906e-3,
                        metallic_radius_m: Some(146e-12),
                    },
                ),
                (
                    "Mo",
                    A {
                        mass_kg_per_mol: 95.95e-3,
                        metallic_radius_m: Some(139e-12),
                    },
                ),
                (
                    "Ag",
                    A {
                        mass_kg_per_mol: 107.87e-3,
                        metallic_radius_m: Some(144e-12),
                    },
                ),
                (
                    "Sn",
                    A {
                        mass_kg_per_mol: 118.71e-3,
                        metallic_radius_m: Some(151e-12),
                    },
                ),
                (
                    "Hf",
                    A {
                        mass_kg_per_mol: 178.49e-3,
                        metallic_radius_m: Some(159e-12),
                    },
                ),
                (
                    "Ta",
                    A {
                        mass_kg_per_mol: 180.95e-3,
                        metallic_radius_m: Some(146e-12),
                    },
                ),
                (
                    "W",
                    A {
                        mass_kg_per_mol: 183.84e-3,
                        metallic_radius_m: Some(139e-12),
                    },
                ),
                (
                    "Pb",
                    A {
                        mass_kg_per_mol: 207.2e-3,
                        metallic_radius_m: Some(175e-12),
                    },
                ),
            ];
            for (s, a) in rows {
                t.insert(*s, *a);
            }
            t
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn iron_mass_is_well_known() {
        let m = PeriodicTable::atomic_mass("Fe").unwrap();
        assert_relative_eq!(m, 55.845e-3, max_relative = 1e-4);
    }

    #[test]
    fn unknown_element_returns_none() {
        assert!(PeriodicTable::atomic_mass("Xx").is_none());
    }

    #[test]
    fn noble_nonmetal_radius_is_none() {
        assert!(PeriodicTable::atomic_radius("O").is_none());
    }

    #[test]
    fn iron_has_radius() {
        let r = PeriodicTable::atomic_radius("Fe").unwrap();
        assert_relative_eq!(r, 126e-12, max_relative = 1e-3);
    }
}
