//! Example: rubber elasticity for a crosslinked elastomer.
//!
//! Computes:
//! 1. The Arruda-Boyce 8-chain Cauchy stress as a function of stretch
//!    for `N = 8` segments and a network shear modulus `G = 0.4 MPa`,
//!    demonstrating the steep upturn near the finite-extensibility
//!    limit (`λ → √N`).
//! 2. The neo-Hookean limit at small stretches.
//! 3. A Worm-Like-Chain (Marko-Siggia) force-extension curve for a
//!    semi-flexible biopolymer.
//!
//! Run with `cargo run -p rubber-elasticity`.

use tpt_mat_polymer::{
    arruda_boyce_true_stress, stress_strain, wlc_force_extension, ChainModel, PolymerModel,
};

fn main() {
    println!("Arruda-Boyce 8-chain rubber elasticity");
    println!("======================================");
    let n = 8_u32;
    let g_pa = 0.4e6_f64;
    println!("  N = {n} segments, G = {g_pa:.1e} Pa\n");

    println!("  λ         σ_AB (Pa)       σ_neoHookean (Pa)");
    println!("  -------   ------------    -----------------");
    let model = PolymerModel {
        chain: ChainModel::ArrudaBoyce {
            n_segments: n,
            shear_modulus: g_pa,
        },
        crosslink_density: 0.0,
    };
    for &lambda in &[1.01, 1.1, 1.5, 2.0, 2.5, 2.7, 2.8] {
        let s_ab = stress_strain(&model, lambda);
        let s_nh = 2.0 * g_pa * (lambda * lambda - 1.0 / (lambda * lambda * lambda));
        println!("  {lambda:>5.2}     {s_ab:>12.3e}    {s_nh:>12.3e}");
    }

    println!(
        "\n  At λ → √N = {:.4} (finite extensibility),",
        (n as f64).sqrt()
    );
    println!("  Arruda-Boyce stress diverges via the inverse-Langevin;");
    println!("  neo-Hookean stays finite (visible above as λ → 2.83).");

    println!("\nWorm-Like-Chain (Marko-Siggia) force-extension");
    println!("==============================================");
    let l_p = 1.5e-9_f64; // 1.5 nm persistence (dsDNA-like)
    let l_c = 50.0e-9_f64; // 50 nm contour
    println!("  l_p = {l_p:.1e} m, L_c = {l_c:.1e} m\n");
    println!("  x/L_c     F (1/k_B T · 1/m, scaled)");
    println!("  -------   --------------------------");
    for &frac in &[0.1, 0.3, 0.5, 0.7, 0.9, 0.97, 0.99] {
        let stretch = frac * l_c;
        let f = wlc_force_extension(l_p, l_c, stretch);
        println!("  {frac:>5.2}     {f:>24.3e}");
    }

    // Demonstrate the public arruda_boyce_true_stress API directly.
    println!("\nDirect arruda_boyce_true_stress(N=4, G=1 MPa, λ=3.0):");
    let s = arruda_boyce_true_stress(4, 1.0e6, 3.0);
    println!("  σ = {s:.3e} Pa");
}
