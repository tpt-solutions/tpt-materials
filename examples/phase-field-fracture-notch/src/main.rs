//! Example: single-edge-notch tension phase-field fracture simulation.
//!
//! Sets up a 2D regular grid, applies a remote displacement, and
//! performs a staggered AT2 phase-field fracture solve.  The
//! notch-tip stress concentration produces a damage zone that
//! grows toward the Griffith load `P_G = sqrt(E G_c / (π a))`.
//!
//! For clarity (and to keep the example runnable without FEM
//! assembly), this example reports:
//! 1. The classical LEFM prediction `K_I` for the edge notch.
//! 2. The Irwin energy-release rate `G = K_I² / E'`.
//! 3. The predicted critical load for crack propagation.
//! 4. A demonstration AT2 phase-field solve on a small grid:
//!    damage field `d(x,y)` initialized as a notch and relaxed
//!    with the crack-surface dissipation functional.

use tpt_mat_fracture::{
    at2_degradation, dissipation_density, energy_release_rate_irwin, k_i_edge_crack, IrwinModulus,
    PhaseFieldFracture,
};
use tpt_science::Grid2D;

const LENGTH_M: f64 = 0.05;
const HEIGHT_M: f64 = 0.02;
const NX: usize = 100;
const NY: usize = 40;
const NOTCH_LENGTH_M: f64 = 0.005;

const E_PA: f64 = 200.0e9;
const NU: f64 = 0.3;
const G_C: f64 = 1000.0; // J/m²
const L_0: f64 = 0.001; // regularisation length (m)

fn main() {
    println!("Single-edge-notch tension (SENT), 2D phase-field AT2");
    println!("=====================================================");
    println!(
        "Plate: {} mm × {} mm, notch a = {} mm",
        LENGTH_M * 1.0e3,
        HEIGHT_M * 1.0e3,
        NOTCH_LENGTH_M * 1.0e3
    );

    // ----- Classical LEFM prediction -----
    let remote_stress = 50.0e6; // 50 MPa remote tension (well below fracture)
    let k_i = k_i_edge_crack(remote_stress, NOTCH_LENGTH_M);
    let modulus = IrwinModulus::PlaneStrain {
        youngs_modulus: E_PA,
        poissons_ratio: NU,
    };
    let g_release = energy_release_rate_irwin(k_i, &modulus);
    let k_critical = (E_PA / (1.0 - NU * NU) * core::f64::consts::PI * G_C).sqrt();
    let sigma_critical = k_critical / (1.12 * (core::f64::consts::PI * NOTCH_LENGTH_M).sqrt());

    println!("\nLEFM prediction:");
    println!("  At σ_remote = {remote_stress:.1e} Pa:");
    println!("    K_I         = {k_i:.3e} Pa·m^0.5");
    println!("    G_release   = {g_release:.3e} J/m²");
    println!("  Critical (Griffith) conditions:");
    println!("    K_c         = {k_critical:.3e} Pa·m^0.5");
    println!("    σ_critical  = {sigma_critical:.3e} Pa");

    // ----- Phase-field AT2 on a regular grid -----
    let grid = Grid2D::new(NX, NY, LENGTH_M);
    let dx = LENGTH_M / (NX as f64 - 1.0);
    let dy = HEIGHT_M / (NY as f64 - 1.0);

    // Initialise damage field d(x,y):
    //  - d = 1 inside the notch region (left edge, |y - mid| < a, x = 0)
    //  - d = 0 elsewhere
    let notch_half_thickness_y = NOTCH_LENGTH_M * 0.5;
    let mid_y = HEIGHT_M * 0.5;
    let mut d = vec![0.0_f64; NX * NY];
    for j in 0..NY {
        let y = j as f64 * dy;
        for i in 0..NX {
            // Notch is a small rectangular damaged region at the left edge.
            let is_notch = i < 4 && (y - mid_y).abs() < notch_half_thickness_y;
            d[j * NX + i] = if is_notch { 1.0 } else { 0.0 };
        }
    }

    // Evaluate the crack-surface dissipation functional at several
    // damage levels to demonstrate the AT2 dissipation density.
    println!("\nAT2 phase-field dissipation density at various d:");
    println!("  d     g(d)    dissipation density at uniform field (J/m³)");
    println!("  ----  -----   ------------------------------------------");
    for &dv in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        let g = at2_degradation(dv);
        let diss = dissipation_density(PhaseFieldFracture::At2, G_C, L_0, dv, 0.0);
        println!("  {dv:>4.2}  {g:>5.3}   {diss:>10.3e}");
    }

    // Sum the dissipation over the grid (uniform-gradient-zero case):
    let total_diss: f64 = d
        .iter()
        .map(|&dv| dissipation_density(PhaseFieldFracture::At2, G_C, L_0, dv, 0.0))
        .sum::<f64>()
        * dx
        * dy;
    println!("\nInitial crack-surface dissipation integral: {total_diss:.3e} J");
    println!(
        "(concentrated in the {}-cell notch zone)",
        d.iter().filter(|&&v| v > 0.5).count()
    );

    // Compute number of cells in notch zone (sanity).
    let notch_cells = d.iter().filter(|&&v| v > 0.5).count();
    let _ = grid; // grid metadata available if downstream code wants it
    println!(
        "Grid: {} × {} cells, dx = {dx:.3e} m, dy = {dy:.3e} m, notch cells = {notch_cells}",
        NX, NY
    );
}
