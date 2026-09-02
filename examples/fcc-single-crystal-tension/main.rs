//! Phase 1 data-only scaffold for the FCC single-crystal tension example.
//!
//! The full solver (assemble stiffness → Newton-Raphson → update
//! stress/hardening/lattice rotation) lands in Phase 2 with
//! `tpt-mat-crystal-plasticity` and `tpt-fem`.
//!
//! This file just prints the orientation and the 12 FCC slip systems so
//! `cargo run --example fcc_single_crystal_tension` produces visible
//! output today.

use tpt_mat_constants::PhysicalConstants;
use tpt_mat_core::{CrystalOrientation, OrientationRepresentation};
use tpt_mat_crystallography::{CrystalStructure, SlipSystem};

fn main() {
    println!("== fcc single-crystal tension (Phase 1 scaffold) ==");
    println!(
        "gas constant  R = {:.6e} J/(mol·K)",
        PhysicalConstants::GAS_CONSTANT
    );
    println!(
        "Avogadro      N_A = {:.6e} /mol",
        PhysicalConstants::AVOGADRO
    );

    let fcc = CrystalStructure::FCC;
    let slips = fcc.slip_systems();
    println!("\nFCC slip systems ({} total):", slips.len());
    for (i, s) in slips.iter().enumerate() {
        println!(
            "  {:>2}. plane {:?}  direction {:?}  CRSS = {:.1} MPa",
            i + 1,
            s.plane_normal,
            s.slip_direction,
            s.critical_resolved_shear_stress
        );
    }

    // Identity orientation: tensile axis = sample X.
    let g = CrystalOrientation::identity();
    let e = g.to_representation(OrientationRepresentation::EulerBunge);
    println!(
        "\nIdentity orientation Euler-Bunge angles: phi1={:.4} Phi={:.4} phi2={:.4}",
        e[0], e[1], e[2]
    );

    // Sample-level tensile axis.
    let tensile_axis = [1.0_f64, 0.0, 0.0];
    let m = SlipSystem::max_schmid_factor(&slips, tensile_axis);
    println!("\nMax Schmid factor for tension along {tensile_axis:?} = {m:.4}");
}
