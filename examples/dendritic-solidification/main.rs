//! Dendritic solidification example (Phase 3).
//!
//! Seeds a small solid nucleus in a supercooled liquid and runs the
//! Kobayashi phase-field solver for a fixed number of steps.  Prints
//! the solid fraction, tip velocity, and arm-spacing estimate at the
//! end of the run.

use tpt_mat_phase_field::BulkEnergy;
use tpt_mat_solidification::{AnisotropyModel, SolidificationSolver};
use tpt_science::Grid2D;

fn main() {
    println!("== dendritic solidification (Phase 3) ==");
    let n = 64;
    let grid = Grid2D::new(n, n, 1.0);

    // Initial order parameter: a 4-cell-radius solid seed in the centre.
    let mut eta = vec![0.0_f64; n * n];
    let cx = n as f64 / 2.0;
    let cy = n as f64 / 2.0;
    let seed_radius = 4.0_f64;
    for i in 0..n {
        for j in 0..n {
            let d = ((i as f64 - cy).powi(2) + (j as f64 - cx).powi(2)).sqrt();
            eta[i * n + j] = if d < seed_radius { 1.0 } else { 0.0 };
        }
    }
    let solid0: f64 = eta.iter().filter(|&&v| v > 0.5).count() as f64;
    println!("initial solid cells: {solid0}");

    let anisotropy = AnisotropyModel {
        strength: 0.04,
        mode: tpt_mat_solidification::AnisotropyMode::Cubic4Fold,
    };
    let mut solver = SolidificationSolver::new(
        grid,
        &eta,
        anisotropy,
        BulkEnergy::DoubleWell { well_depth: 1.0 },
        0.005,
        1.0,
        1.0,
    )
    .unwrap();

    let snap0 = solver.snapshot();
    println!("t = 0: solid_fraction = {:.4}", snap0.solid_fraction);

    let n_steps = 100;
    for step in 0..n_steps {
        solver.step().unwrap();
        if step % 20 == 19 {
            let s = solver.snapshot();
            println!(
                "after step {step}: solid_fraction = {:.4} (ΔE grown)",
                s.solid_fraction
            );
        }
    }
    let snap = solver.snapshot();
    println!("---");
    println!("final solid_fraction = {:.4}", snap.solid_fraction);
    println!("tip velocity (cells / time) ≈ {:.4}", snap.tip_velocity);
    println!(
        "primary arm spacing (cells) = {:.2}",
        snap.primary_arm_spacing
    );
    println!(
        "secondary arm spacing (cells) = {:.2}",
        snap.secondary_arm_spacing
    );

    // Note: tip velocity estimate uses the full domain length as a
    // length scale, which is a coarse proxy for true tip speed in
    // a periodic domain.
}
