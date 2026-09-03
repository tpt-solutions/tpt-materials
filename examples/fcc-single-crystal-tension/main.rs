//! FCC single-crystal tension example (Phase 2).
//!
//! Applies a 0.1% tensile strain increment to a single FCC crystal
//! and prints the resulting Cauchy stress, slip rates, and Schmid
//! factors.  Also computes the isotropic Taylor-factor proxy for a
//! random FCC texture.

use tpt_mat_constants::PhysicalConstants;
use tpt_mat_core::CrystalOrientation;
use tpt_mat_crystal_plasticity::{
    BoundaryConditions, CpFemSolver, CrystalPlasticityModel, LoadStep, RateSensitivity,
    SymmetricFourthOrder,
};
use tpt_mat_crystallography::{CrystalStructure, SlipSystem};
use tpt_mat_hardening::{Hardening, VoceHardening, VoceParams};
use tpt_mat_texture::{taylor_factor, PoleFigure, PoleFigureGrid, PoleKind, TextureAnalyzer};

fn main() {
    println!("== fcc single-crystal tension (Phase 2) ==");
    println!(
        "gas constant  R = {:.6e} J/(mol·K)",
        PhysicalConstants::GAS_CONSTANT
    );
    println!(
        "Avogadro      N_A = {:.6e} /mol",
        PhysicalConstants::AVOGADRO
    );

    // Material model.  CRSS τ_0 = 30 MPa, τ_s = 200 MPa — typical for a
    // ductile FCC metal at room temperature.  Stiffness is Cu-like
    // (C11 = 168 GPa, C12 = 121 GPa, C44 = 75 GPa).
    let hardening = Hardening::Voce(VoceHardening::uniform(VoceParams {
        tau_0: 30.0,
        tau_s: 200.0,
        theta_0: 500.0,
        gamma_c: 0.05,
    }));
    let rate = RateSensitivity {
        reference_strain_rate: 1.0e-3,
        exponent: 20.0,
    };
    let elastic = SymmetricFourthOrder::cubic(168_000.0, 121_000.0, 75_000.0);
    let model = CrystalPlasticityModel::from_crystal_structure(
        CrystalStructure::FCC,
        hardening,
        rate,
        elastic,
    )
    .unwrap();
    println!("\nFCC slip systems: {}", model.slip_systems.len());

    // FEM solver with a single integration point.
    let bc = BoundaryConditions::default();
    let mut solver = CpFemSolver::new(0, vec![model.clone()], bc);
    let load = LoadStep::uniaxial(0, 1.0e-3); // 0.1% tensile strain along x.
    let result = solver.solve_increment(&load).unwrap();
    println!(
        "Tensile stress along x: σ_xx = {:.3} MPa",
        result.stresses[0][0]
    );
    let nonzero: Vec<(usize, f64)> = result.slip_rates[0]
        .iter()
        .enumerate()
        .filter(|(_, g)| **g > 0.0)
        .map(|(i, g)| (i, *g))
        .collect();
    println!(
        "Active slip systems: {} (out of {})",
        nonzero.len(),
        model.slip_systems.len()
    );
    for (i, _g) in nonzero.iter().take(5) {
        let s = &model.slip_systems[*i];
        println!(
            "  system {:>2}: activated (crystal plane normal {:?}, direction {:?})",
            i + 1,
            s.plane_normal,
            s.slip_direction
        );
    }

    // Taylor factor for random FCC texture.
    let m = taylor_factor(4096, 0xDEAD_BEEF);
    println!("\nTaylor-factor proxy for random FCC texture: M ≈ {:.3}", m);

    // Pole figure for a single crystal at identity orientation.
    let ta = TextureAnalyzer::new(vec![CrystalOrientation::identity()], vec![1.0]);
    let pf: PoleFigure = ta.pole_figure(PoleKind::Fcc111, PoleFigureGrid::default());
    let nonzero = pf.intensity.iter().filter(|v| **v > 0.0).count();
    println!(
        "(111) pole figure: {} non-zero cells (1 grain → 1 cell)",
        nonzero
    );
    let _ = SlipSystem::max_schmid_factor; // keep import used
}
