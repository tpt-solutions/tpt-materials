//! Example: Eshelby (1957) equivalent-inclusion for a stiff
//! spherical particle in an aluminium matrix.
//!
//! Demonstrates:
//!
//! - [`tpt_mat_homogenization::eshelby_spherical`] — Eshelby tensor.
//! - [`tpt_mat_homogenization::dilute_strain_concentration`] — the
//!   dilute strain-concentration tensor.
//! - [`tpt_mat_composite_micro::mori_tanaka`] — the Mori–Tanaka
//!   mean-field homogenization.
//!
//! Run with: `cargo run --release -p eshelby-inclusion`.

use tpt_mat_composite_micro::{dilute_estimate, eshelby_for_isotropic_sphere, mori_tanaka};
use tpt_mat_crystal_plasticity::SymmetricFourthOrder;
use tpt_mat_homogenization::dilute_strain_concentration;
use tpt_math_linalg_fixed::Vec6;

fn main() {
    // Aluminium matrix: E = 70 GPa, ν = 0.33.
    let e_m = 70_000.0_f64;
    let nu_m = 0.33_f64;
    // SiC-like inclusion: E = 400 GPa, ν = 0.25.
    let e_i = 400_000.0_f64;
    let nu_i = 0.25_f64;
    let c_matrix = SymmetricFourthOrder::isotropic(e_m, nu_m);
    let c_inclusion = SymmetricFourthOrder::isotropic(e_i, nu_i);

    // Eshelby tensor for a spherical inclusion.
    let s = eshelby_for_isotropic_sphere(e_m, nu_m);
    println!("Eshelby tensor (spherical, ν = {nu_m}):");
    println!("  S_h = {:.4} (hydrostatic)", s.s_hydro);
    println!("  S_d = {:.4} (deviatoric)", s.s_dev);

    // Dilute strain-concentration tensor.
    let a = dilute_strain_concentration(&c_matrix, &c_inclusion, s);
    println!();
    println!("Dilute strain-concentration tensor A:");
    for i in 0..6 {
        let row: Vec<String> = (0..6).map(|j| format!("{:>8.4}", a.data[i][j])).collect();
        println!("  [{}]", row.join(", "));
    }

    // Apply far-field strain ε_xx = 1.0, others 0.
    let eps_inf = Vec6::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let eps_inc = a.apply(eps_inf);
    println!();
    println!("For ε^∞ = (1, 0, 0, 0, 0, 0), inclusion strain ε^inc:");
    println!(
        "  ε_xx = {:.4} (< 1.0 because inclusion is stiffer)",
        eps_inc[0]
    );
    println!("  ε_yy = ε_zz = {:.4} (lateral contraction)", eps_inc[1]);

    // Mori–Tanaka: sweep inclusion volume fraction.
    println!();
    println!("Mori–Tanaka effective stiffness (C_eff[0][0]) as a function of f_inclusion:");
    println!("  f_inc    C_eff[0][0]    C_eff[1][1]    C_eff[3][3]");
    println!("  -----    ----------    ----------    ----------");
    for i in 0..=10 {
        let f = i as f64 / 10.0;
        let mt = mori_tanaka(&c_matrix, &c_inclusion, f, s);
        println!(
            "  {:>5.2}    {:>10.0}    {:>10.0}    {:>10.0}",
            f, mt.c_eff.data[0][0], mt.c_eff.data[1][1], mt.c_eff.data[3][3]
        );
    }

    // Dilute estimate at f = 0.1.
    let f = 0.1_f64;
    let c_dilute = dilute_estimate(&c_matrix, &c_inclusion, f, s);
    println!();
    println!(
        "Dilute estimate at f = {f}: C_eff[0][0] = {:.0} MPa",
        c_dilute.data[0][0]
    );
    println!(
        "Mori–Tanaka at f = {f}: C_eff[0][0] = {:.0} MPa",
        mori_tanaka(&c_matrix, &c_inclusion, f, s).c_eff.data[0][0]
    );
}
