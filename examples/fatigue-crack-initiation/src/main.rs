//! Example: microstructural-fatigue crack-initiation prediction.
//!
//! Builds a small 5-grain polycrystal RVE, runs a cycle-by-cycle
//! Findley FIP scan over the (synthetic) accumulated-shear history
//! and identifies the most likely crack-initiation site.
//!
//! Run with: `cargo run --release -p fatigue-crack-initiation`.

use tpt_mat_crystal_plasticity::SymmetricFourthOrder;
use tpt_mat_fatigue_micro::{
    predict_crack_initiation, FatigueCriterion, FatigueCriterion::Findley,
};
use tpt_mat_rve::Rve;
use tpt_math_linalg_fixed::Vec6;

fn main() {
    let n_grains = 5;
    // Build a 5-grain RVE with random orientations.
    let mut rve = Rve::empty();
    for i in 0..n_grains {
        let angle = (i as f64) * core::f64::consts::PI / 7.0;
        let c = angle.cos();
        let s = angle.sin();
        let orientation = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let c_stiff = SymmetricFourthOrder::cubic(200_000.0, 100_000.0, 50_000.0);
        rve.grains.push(tpt_mat_rve::RveGrain::new(
            format!("g{i}"),
            1.0 / n_grains as f64,
            orientation,
            c_stiff,
        ));
    }
    // Synthetic accumulated-shear history — most grains carry a small
    // amount of slip, grain 2 carries ~10× more.
    let mut gamma_acc = vec![vec![0.001_f64; 12]; n_grains];
    gamma_acc[2] = vec![0.05_f64; 12];

    let criterion: FatigueCriterion = Findley { k: 0.3 };
    let res = predict_crack_initiation(
        &gamma_acc, criterion, 1000.0, // cycles in the supplied history
        0.5,    // Δγ per cycle (Coffin–Manson driver)
    );

    println!("Microstructural fatigue crack-initiation prediction");
    println!("==================================================");
    println!("RVE: {n_grains} grains, cubic-stiffness isotropic");
    println!("Criterion: Findley (k = 0.3)");
    println!();
    println!("FIP field:");
    for (i, &fip) in res.fip_field.iter().enumerate() {
        println!("  grain {i}: FIP = {fip:.3}");
    }
    println!();
    println!("Critical grain:        {}", res.critical_grain);
    println!("Critical FIP:          {:.3}", res.critical_fip);
    println!("Cycles-to-initiation:  {:.0}", res.cycles_to_initiation);
}
