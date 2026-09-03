//! Spinodal decomposition example (Phase 3).
//!
//! Initialises a 64×64 concentration field with a small sinusoidal
//! perturbation around the spinodal composition, advances the
//! Cahn–Hilliard equation for a few hundred steps, and reports the
//! final free energy / interface area.  A free-energy decrease over
//! time is the canonical verification that the discretisation is
//! dissipative.

use tpt_mat_phase_field::{
    BulkEnergy, PhaseFieldModel, PhaseFieldSolver, RegularSolutionParams,
};
use tpt_science::Grid2D;

fn main() {
    println!("== spinodal decomposition (Phase 3) ==");
    let n = 64;
    let grid = Grid2D::new(n, n, 1.0);
    let mut c = vec![0.5_f64; n * n];
    // Initial concentration: a 5% perturbation around the spinodal
    // composition c = 0.5.  For omega = 4, RT = 1 the spinodal
    // condition `f''(c) = -2Ω + RT / (c(1-c)) < 0` is satisfied
    // around c = 0.5, so the perturbation grows (decomposition).
    for (idx, v) in c.iter_mut().enumerate() {
        let i = idx / n;
        let j = idx % n;
        let phi = 2.0 * std::f64::consts::PI
            * ((i as f64) * 0.05 + (j as f64) * 0.07);
        *v = 0.5 + 0.05 * phi.sin() + 0.03 * (2.0 * phi).cos();
    }
    let mean: f64 = c.iter().sum::<f64>() / c.len() as f64;
    println!("initial mean concentration = {mean:.6}");

    let mut solver = PhaseFieldSolver::new_cahn_hilliard(
        grid,
        0.001,
        1.0,
        1.0,
        BulkEnergy::RegularSolution(RegularSolutionParams { omega: 4.0, rt: 1.0 }),
        &c,
    )
    .unwrap();

    let e0 = solver.snapshot().free_energy;
    let a0 = solver.snapshot().interface_area;
    let n_steps = 300;
    solver.step_many(n_steps).unwrap();
    let snap = solver.snapshot();
    println!(
        "after {n_steps} steps: free energy {e0:.3} → {:.3} (Δ = {:.3e})",
        snap.free_energy,
        snap.free_energy - e0
    );
    println!(
        "interface area: {a0:.1} → {:.1}",
        snap.interface_area
    );

    let mean_final: f64 = solver
        .concentration
        .iter()
        .sum::<f64>()
        / solver.concentration.len() as f64;
    println!("final mean concentration = {mean_final:.6}");

    assert!(snap.free_energy <= e0 + 1e-6, "free energy must not increase");
    // Neumann-flux boundary scheme conserves mean composition to round-off
    // times grid size; allow a small numerical drift.
    assert!(
        (mean_final - mean).abs() < 1e-3,
        "Cahn-Hilliard must conserve mean composition (drift = {:.2e})",
        (mean_final - mean).abs()
    );
    println!("OK: free energy non-increasing and mean concentration conserved.");
}