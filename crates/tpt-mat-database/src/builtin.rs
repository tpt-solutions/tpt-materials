//! Bundled MIT-clean materials data set.

use super::provenance::DataSource;
use super::record::{
    Composition, Electrical, Mechanical, MaterialRecord, Thermal,
};

/// Build the bundled material records.
///
/// All values are textbook reference numbers or values that the
/// project can publish under MIT/Apache-2.0.  Each record carries
/// the [`DataSource`] provenance required by spec §9.
pub fn builtin_records() -> Vec<MaterialRecord> {
    vec![
        MaterialRecord {
            name: "AISI 304 (stainless steel)".to_string(),
            composition: Composition {
                fe_balance: 0.70,
                cr: 0.19,
                ni: 0.10,
                c: 0.0008,
                ..Default::default()
            },
            mechanical: Mechanical {
                youngs_modulus_gpa: 200.0,
                yield_strength_mpa: 215.0,
                ultimate_tensile_strength_mpa: 505.0,
                density_kg_per_m3: 8000.0,
                poissons_ratio: 0.29,
            },
            thermal: Thermal {
                thermal_conductivity_w_per_mk: 16.2,
                specific_heat_j_per_kgk: 500.0,
                coefficient_thermal_expansion_per_k: 17.3e-6,
                melting_point_k: 1670.0,
            },
            electrical: Electrical {
                electrical_resistivity_ohm_m: 7.2e-7,
            },
            sources: vec![DataSource::Textbook(
                "ASM Handbook Vol. 1 (1990)".to_string(),
            )],
        },
        MaterialRecord {
            name: "AISI 4140 (Cr-Mo steel)".to_string(),
            composition: Composition {
                fe_balance: 0.97,
                cr: 0.01,
                mo: 0.002,
                c: 0.004,
                ..Default::default()
            },
            mechanical: Mechanical {
                youngs_modulus_gpa: 205.0,
                yield_strength_mpa: 655.0,
                ultimate_tensile_strength_mpa: 1020.0,
                density_kg_per_m3: 7850.0,
                poissons_ratio: 0.30,
            },
            thermal: Thermal {
                thermal_conductivity_w_per_mk: 42.6,
                specific_heat_j_per_kgk: 473.0,
                coefficient_thermal_expansion_per_k: 12.3e-6,
                melting_point_k: 1570.0,
            },
            electrical: Electrical {
                electrical_resistivity_ohm_m: 2.2e-7,
            },
            sources: vec![DataSource::Textbook(
                "ASM Handbook Vol. 1 (1990)".to_string(),
            )],
        },
        MaterialRecord {
            name: "Al 6061-T6".to_string(),
            composition: Composition {
                al_balance: 0.96,
                mg: 0.01,
                si: 0.006,
                ..Default::default()
            },
            mechanical: Mechanical {
                youngs_modulus_gpa: 69.0,
                yield_strength_mpa: 276.0,
                ultimate_tensile_strength_mpa: 310.0,
                density_kg_per_m3: 2700.0,
                poissons_ratio: 0.33,
            },
            thermal: Thermal {
                thermal_conductivity_w_per_mk: 167.0,
                specific_heat_j_per_kgk: 896.0,
                coefficient_thermal_expansion_per_k: 23.0e-6,
                melting_point_k: 855.0,
            },
            electrical: Electrical {
                electrical_resistivity_ohm_m: 4.0e-8,
            },
            sources: vec![DataSource::Textbook(
                "ASM Handbook Vol. 2 (1990)".to_string(),
            )],
        },
        MaterialRecord {
            name: "Ti-6Al-4V".to_string(),
            composition: Composition {
                ti_balance: 0.90,
                al: 0.06,
                v: 0.04,
                ..Default::default()
            },
            mechanical: Mechanical {
                youngs_modulus_gpa: 113.0,
                yield_strength_mpa: 880.0,
                ultimate_tensile_strength_mpa: 950.0,
                density_kg_per_m3: 4430.0,
                poissons_ratio: 0.34,
            },
            thermal: Thermal {
                thermal_conductivity_w_per_mk: 6.7,
                specific_heat_j_per_kgk: 526.0,
                coefficient_thermal_expansion_per_k: 8.6e-6,
                melting_point_k: 1878.0,
            },
            electrical: Electrical {
                electrical_resistivity_ohm_m: 1.7e-6,
            },
            sources: vec![DataSource::Textbook(
                "ASM Handbook Vol. 2 (1990)".to_string(),
            )],
        },
        MaterialRecord {
            name: "NMC811 (cathode)".to_string(),
            composition: Composition {
                o_balance: 0.4,
                ni: 0.8,
                mn: 0.1,
                co: 0.1,
                li: 1.0,
                ..Default::default()
            },
            mechanical: Mechanical {
                youngs_modulus_gpa: 175.0,
                yield_strength_mpa: 250.0,
                ultimate_tensile_strength_mpa: 350.0,
                density_kg_per_m3: 4770.0,
                poissons_ratio: 0.25,
            },
            thermal: Thermal {
                thermal_conductivity_w_per_mk: 3.5,
                specific_heat_j_per_kgk: 900.0,
                coefficient_thermal_expansion_per_k: 12.0e-6,
                melting_point_k: 1100.0,
            },
            electrical: Electrical {
                electrical_resistivity_ohm_m: 1.0e-3,
            },
            sources: vec![DataSource::Textbook(
                "Daniel et al., J. Power Sources 2014".to_string(),
            )],
        },
    ]
}