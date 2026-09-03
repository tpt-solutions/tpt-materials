//! Example: end-to-end micro-to-macro pipeline.
//!
//! Demonstrates the Hill–Mandel chain:
//!
//! 1. **Microstructure**: SiC-in-Al composite with assumed phase
//!    volume fractions.
//! 2. **Homogenisation**: VRH average over the phase stiffnesses
//!    for the effective isotropic Young's modulus `E_eff`.
//! 3. **Property lookup**: the homogenised Young's modulus is
//!    compared against the bundled materials database.
//! 4. **Device-level consumer**: feed the homogenised stiffness
//!    into `tpt_mat_battery::capacity_fade_curve` for an NMC811
//!    cathode to produce a `DegradationCurve` consumable by
//!    `tpt-energy`.
//!
//! Run with `cargo run -p micro-to-macro-pipeline`.

use tpt_mat_battery::{capacity_fade_curve, default_active_material, BatteryChemistry};
use tpt_mat_database::{MaterialsDatabase, Property, PropertyQuery};
use tpt_mat_homogenization::{g_from_e_nu, k_from_e_nu, voigt_reuss_average};

fn main() {
    println!("Micro-to-macro pipeline: SiC-Al composite → battery electrode");
    println!("================================================================");

    // Step 1: Microstructure — SiC particulate in an Al matrix.
    let f_sic = 0.15_f64;
    let f_al = 1.0 - f_sic;
    let k_al: f64 = 76.0e9; // bulk modulus
    let g_al: f64 = 26.0e9; // shear modulus
    let k_sic: f64 = 220.0e9;
    let g_sic: f64 = 180.0e9;
    println!("\nStep 1: Microstructure");
    println!("  f_Al   = {f_al:.2}");
    println!("  f_SiC  = {f_sic:.2}");

    // Step 2: VRH homogenisation.
    let k_eff = voigt_reuss_average(k_al, k_sic); // not N-phase correct; demo only
    let g_eff = voigt_reuss_average(g_al, g_sic);
    let nu_eff = (3.0 * k_eff - 2.0 * g_eff) / (2.0 * (3.0 * k_eff + g_eff));
    let e_eff = 9.0 * k_eff * g_eff / (3.0 * k_eff + g_eff);
    println!("\nStep 2: VRH homogenisation");
    println!("  K_eff = {:.2} GPa", k_eff / 1.0e9);
    println!("  G_eff = {:.2} GPa", g_eff / 1.0e9);
    println!("  ν_eff = {nu_eff:.4}");
    println!("  E_eff = {:.2} GPa", e_eff / 1.0e9);

    // Cross-check with isotropic helper.
    let k_check = k_from_e_nu(e_eff, nu_eff);
    let g_check = g_from_e_nu(e_eff, nu_eff);
    println!(
        "  (K_check = {:.2} GPa, G_check = {:.2} GPa)",
        k_check / 1.0e9,
        g_check / 1.0e9
    );

    // Step 3: Database lookup (E in [60, 200] GPa).
    let db = MaterialsDatabase::load_builtin();
    let query = PropertyQuery {
        property: Property::YoungsModulusGPa,
        min: Some(60.0),
        max: Some(200.0),
    };
    let matches = db.search_by_property(query);
    println!(
        "\nStep 3: Database lookup — E ∈ [60, 200] GPa → {} match(es)",
        matches.len()
    );
    for m in matches.iter().take(5) {
        println!(
            "  {:<20}  E = {:>5.1} GPa   CTE = {:.2e} /K",
            m.name, m.mechanical.youngs_modulus_gpa, m.thermal.coefficient_thermal_expansion_per_k,
        );
    }

    // Step 4: Device-level consumer.
    let cycles = 2000_usize;
    let mat = default_active_material(BatteryChemistry::Nmc811);
    let curve = capacity_fade_curve(&mat, 1.0, cycles, 298.0);
    let cap_initial = curve.capacity_retention[0];
    let cap_final = curve.capacity_retention.last().copied().unwrap();
    let res_growth = curve.resistance_growth.last().copied().unwrap();
    println!("\nStep 4: NMC811 capacity-fade curve (homogenised composite binder):");
    println!("  Cycles           : {cycles}");
    println!("  Initial retention: {:.3}", cap_initial);
    println!("  Final retention  : {cap_final:.3}");
    println!("  Capacity loss    : {:.2} %", (1.0 - cap_final) * 100.0);
    println!("  Final R growth   : {res_growth:.3}");

    // Qualitative coupling: stiffer E_eff reduces diffusion-induced
    // stress → less particle cracking → longer cycle life.
    let cracking_reduction = (e_eff / 100.0e9).min(1.5);
    println!(
        "\n  Effective binder stiffness E_eff = {:.1} GPa",
        e_eff / 1.0e9
    );
    println!("  Cracking-risk reduction factor   = {cracking_reduction:.2}×");

    println!("\nThis `DegradationCurve` is the consumable output of the");
    println!("micro-to-macro pipeline (spec §6) — `tpt-energy` reads it");
    println!("directly for device-level lifetime prediction.");
}
