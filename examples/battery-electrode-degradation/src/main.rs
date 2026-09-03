//! Example: NMC811 capacity fade over 10 000 cycles.
//!
//! Produces a [`DegradationCurve`] ready to be consumed by the
//! cross-repo `tpt-energy` crate (spec §6).

use tpt_mat_battery::{capacity_fade_curve, default_active_material, BatteryChemistry};

fn main() {
    let material = default_active_material(BatteryChemistry::Nmc811);
    let curve = capacity_fade_curve(&material, 1.0, 10_000, 298.0);

    println!("NMC811 capacity-fade curve, 1 C, 298 K, 10 000 cycles");
    println!("==================================================");
    println!(
        "  cycle {:>5} → retention = {:.3}, resistance growth = {:.2}×",
        curve.cycles[0], curve.capacity_retention[0], curve.resistance_growth[0]
    );
    let mut milestones = vec![10_usize, 100, 500, 1_000, 5_000, 10_000];
    milestones.retain(|n| *n <= curve.cycles.len());
    for n in milestones {
        let i = n - 1;
        println!(
            "  cycle {:>5} → retention = {:.3}, resistance growth = {:.2}×",
            curve.cycles[i], curve.capacity_retention[i], curve.resistance_growth[i]
        );
    }
    let final_cap = *curve.capacity_retention.last().unwrap();
    println!();
    println!("Final capacity retention: {:.1} %", 100.0 * final_cap);
}
