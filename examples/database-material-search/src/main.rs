//! Materials-database property-search example.

use tpt_mat_database::{MaterialsDatabase, Property, PropertyQuery};

fn main() {
    let db = MaterialsDatabase::load_builtin();
    println!(
        "Loaded {} bundled materials; first 3 names:",
        db.len()
    );
    for m in db.materials.iter().take(3) {
        println!("  - {}", m.name);
    }

    let q = PropertyQuery::new(Property::YoungsModulusGPa)
        .min(150.0)
        .max(300.0);
    println!("\nMaterials with Young's modulus in [150, 300] GPa:");
    for h in db.search_by_property(q) {
        println!(
            "  {:<30} E = {:>6.1} GPa   σ_y = {:>6.1} MPa   ρ = {:>5.0} kg/m³",
            h.name,
            h.mechanical.youngs_modulus_gpa,
            h.mechanical.yield_strength_mpa,
            h.mechanical.density_kg_per_m3
        );
    }

    let q = PropertyQuery::new(Property::ThermalConductivityWPerMK).min(40.0);
    println!("\nMaterials with thermal conductivity ≥ 40 W/(m·K):");
    for h in db.search_by_property(q) {
        println!(
            "  {:<30} k = {:>6.2} W/(m·K)   T_m = {:>5.0} K",
            h.name,
            h.thermal.thermal_conductivity_w_per_mk,
            h.thermal.melting_point_k
        );
    }
}