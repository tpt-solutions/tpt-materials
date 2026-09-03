//! Example: LaNi5 metal-hydride PCT isotherms.
//!
//! Computes pressure-composition isotherms at three temperatures
//! and prints the plateau-pressure shift driven by the van 't Hoff
//! relation.  The ΔH_des ≈ 30 kJ/mol H2 used here is the standard
//! textbook value for LaNi5 (endothermic desorption: increasing T
//! raises the plateau pressure).

use tpt_mat_hydrogen_storage::{
    pct_isotherm, van_t_hoff_pressure, HydrideType, HydrogenStorageMaterial, PctParams,
};

fn main() {
    // Standard LaNi5H6 plateau pressure is about 2 bar at 298 K.
    // ΔH_des = +30 kJ/mol (endothermic), ΔS_des ≈ +100 J/(mol·K).
    let params = PctParams {
        p0_at_tref: 2.0e5,
        t_ref_k: 298.0,
        delta_h_j_per_mol: 30_000.0,
        delta_s_j_per_mol_k: 100.0,
        h_to_m_max: 6.0,
        slope_low: 1.0e5,
        slope_up: 1.0e5,
    };

    let lani5 = HydrogenStorageMaterial {
        hydride_type: HydrideType::MetalHydride {
            alloy: "LaNi5".to_string(),
        },
        storage_capacity: 1.4, // wt%
        desorption_enthalpy_j_per_mol: 30_000.0,
        desorption_entropy_j_per_mol_k: 100.0,
    };
    println!(
        "Material: {:?}  capacity = {:.2} wt%",
        lani5.hydride_type, lani5.storage_capacity
    );
    println!();

    println!("Van 't Hoff plateau pressure vs temperature:");
    println!("  T (K)   P_eq (bar)");
    println!("  -----   ----------");
    for &t in &[298.0, 323.0, 348.0, 373.0] {
        let p_pa = van_t_hoff_pressure(&params, t);
        println!("  {t:>5.0}   {:>9.3}", p_pa / 1.0e5);
    }

    println!("\nPCT isotherms (11 points each, H/M ∈ [0, 6]):");
    println!("  H/M     P(298 K)   P(348 K)   P(398 K)");
    println!("        (bar)       (bar)       (bar)");
    println!("  ---     -------    -------    -------");
    let iso_298 = pct_isotherm(&params, 298.0, 11);
    let iso_348 = pct_isotherm(&params, 348.0, 11);
    let iso_398 = pct_isotherm(&params, 398.0, 11);
    for i in 0..iso_298.len() {
        let h = iso_298[i].h_to_m;
        let p1 = iso_298[i].pressure_pa / 1.0e5;
        let p2 = iso_348[i].pressure_pa / 1.0e5;
        let p3 = iso_398[i].pressure_pa / 1.0e5;
        println!("  {h:>4.1}   {p1:>7.3}    {p2:>7.3}    {p3:>7.3}");
    }

    // Show that the plateau widens with temperature (van 't Hoff shift).
    let p_low = van_t_hoff_pressure(&params, 298.0) / 1.0e5;
    let p_high = van_t_hoff_pressure(&params, 398.0) / 1.0e5;
    println!("\nPlateau shift: P_eq(298 K) = {p_low:.3} bar → P_eq(398 K) = {p_high:.3} bar");
    println!("(positive ΔH_des ⇒ dP/dT > 0)");
}
