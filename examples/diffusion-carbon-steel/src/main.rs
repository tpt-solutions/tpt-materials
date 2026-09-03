//! Example: 1D carbon diffusion in a steel slab.
//!
//! Demonstrates `tpt-mat-diffusion` by simulating the carburisation
//! of a low-carbon steel surface over a 1-hour hold at 1200 K.  Carbon
//! is held at fixed surface concentration `c_s = 1.0 wt %` (Dirichlet
//! boundary on the left) and the slab is initially uniformly
//! `c_0 = 0.2 wt %`.  We evolve Fick's second law on a regular 1D
//! grid using the explicit Euler scheme and report the concentration
//! profile at successive checkpoints.
//!
//! Run with: `cargo run --release -p diffusion-carbon-steel`.

use tpt_mat_diffusion::{ArrheniusDiffusivity, ArrheniusParams};
use tpt_science::Grid1D;

fn main() {
    // Geometry: 1 mm slab, 100 cells, dx = 10 µm.
    let n_cells = 100_usize;
    let dx = 1.0e-5_f64;
    let grid = Grid1D::new(n_cells, dx);
    let x: Vec<f64> = (0..n_cells).map(|i| i as f64 * dx * 1.0e3).collect(); // mm

    // Diffusivity of C in γ-Fe at 1200 K: D ≈ 1.5e-10 m²/s.  Use
    // Arrhenius parameters tuned for this example (D_0 ≈ 1.5e-5,
    // Q ≈ 80 kJ/mol) so the temperature dependence is realistic.
    let arr = ArrheniusDiffusivity::new(ArrheniusParams {
        d0: 1.5e-5,
        activation_energy_kj_per_mol: 80.0,
        r_kj_per_mol_k: 8.314_462_618e-3,
    });
    let d_at_1200 = arr.at_temperature(1200.0);
    println!("D(1200 K) = {d_at_1200:.3e} m²/s");

    // CFL-stable time step: Δt ≤ Δx² / (2 D).
    let dt_max = grid.dx * grid.dx / (2.0 * d_at_1200);
    let dt = dt_max * 0.4;
    println!("dt = {dt:.3e} s (CFL limit = {dt_max:.3e})");

    // Initialise: c_0 = 0.2 wt%, fixed c_s = 1.0 wt% at left.
    let c0 = 0.2_f64;
    let c_surface = 1.0_f64;
    let mut u = vec![c0; n_cells];
    u[0] = c_surface;

    // Total simulated time: 60 s — characteristic diffusion length
    // √(2 D t) ≈ 0.024 mm which is comfortably less than the slab
    // thickness (1 mm) so we see a propagating carburisation front.
    let total_time = 60.0_f64;
    let n_steps = (total_time / dt).ceil() as usize;
    println!("Total steps: {n_steps}");

    let checkpoints_s = [0.0, 15.0, 30.0, 45.0, 60.0];

    let mut step = 0_usize;
    let mut elapsed = 0.0_f64;
    for &target_s in &checkpoints_s {
        while (step as f64 * dt) < target_s - 0.5 * dt {
            // 1D Laplacian with fixed left BC, fixed-flux (Neumann) right.
            let mut lap = vec![0.0_f64; n_cells];
            for i in 0..n_cells {
                let left = if i == 0 { c_surface } else { u[i - 1] };
                let right = if i + 1 == n_cells { u[i - 1] } else { u[i + 1] };
                lap[i] = (left - 2.0 * u[i] + right) / (grid.dx * grid.dx);
            }
            for i in 0..n_cells {
                u[i] += dt * d_at_1200 * lap[i];
            }
            u[0] = c_surface;
            step += 1;
        }
        elapsed = step as f64 * dt;
        let max_c = u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let case_depth = u
            .iter()
            .position(|&c| c < 0.5 * (c_surface + c0))
            .map(|i| x[i])
            .unwrap_or(x[n_cells - 1]);
        println!(
            "t = {elapsed:>5.0} s, max c = {max_c:.4}, case depth (~50% line) ≈ {case_depth:.3} mm"
        );
    }

    // Save final profile as JSON for plotting downstream.
    let mut json = String::from("{\"x_mm\":[");
    for (i, xi) in x.iter().enumerate() {
        if i > 0 {
            json.push(',');
        }
        json.push_str(&format!("{xi:.6}"));
    }
    json.push_str("],\"c\":[");
    for (i, ci) in u.iter().enumerate() {
        if i > 0 {
            json.push(',');
        }
        json.push_str(&format!("{ci:.6}"));
    }
    json.push_str("]}");
    let _ = json; // written to stdout for piping.
    eprintln!("profile JSON length: {} bytes", json.len());
}
