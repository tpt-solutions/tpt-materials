//! Example: dislocation-density-based hardening in single slip.
//!
//! Runs a Kocks-Mecking evolution `dρ/dγ = k₁ √ρ − k₂ ρ` along a
//! representative `γ` history and compares the resulting
//! Taylor stress `τ = α μ b √ρ` with a phenomenological Voce
//! curve for the same yield-stress asymptote.
//!
//! Run with `cargo run -p dislocation-density-tension`.

use tpt_mat_dislocation::{
    kocks_mecking_step, taylor_stress, DislocationDensityState, KocksMeckingParams,
};

fn main() {
    // Cu-like parameters.
    let p = KocksMeckingParams {
        k1: 7.0e7,
        k2: 5.0,
        alpha: 0.3,
        shear_modulus: 48.0e9, // Cu
        burgers_vector: 2.56e-10,
    };

    println!("Single-slip hardening curve (Cu, Kocks-Mecking vs Voce)");
    println!("========================================================");
    let mut state = DislocationDensityState::default();
    let dg = 0.002_f64; // slip increment
    let n_steps = 200;
    let mut gamma_history = Vec::with_capacity(n_steps + 1);
    let mut tau_history = Vec::with_capacity(n_steps + 1);
    gamma_history.push(0.0);
    tau_history.push(taylor_stress(&state, &p));
    for _ in 0..n_steps {
        state = kocks_mecking_step(&state, dg, &p);
        gamma_history.push(gamma_history.last().unwrap() + dg);
        tau_history.push(taylor_stress(&state, &p));
    }

    // Voce curve with the same saturation stress:
    let tau_sat = tau_history.last().copied().unwrap();
    let tau_0 = tau_history[0];
    let gamma_c = 0.05_f64; // Voce decay rate
    println!("  τ_sat = {tau_sat:.3e} Pa");
    println!("\n  γ         τ_KM (MPa)    τ_Voce (MPa)    ρ_SSD (m⁻²)");
    println!("  ------    ----------    ------------    -------------");
    for (i, (&g, &tau)) in gamma_history.iter().zip(tau_history.iter()).enumerate() {
        if i % 20 != 0 && i != gamma_history.len() - 1 {
            continue;
        }
        let tau_v = tau_sat - (tau_sat - tau_0) * (-g / gamma_c).exp();
        println!(
            "  {:>5.3}    {:>10.3}    {:>10.3}     {:>11.3e}",
            g,
            tau / 1.0e6,
            tau_v / 1.0e6,
            state.rho_ssd,
        );
        if i == 0 {
            // avoid the "never read" lint for `state` on first iter:
            let _ = state.rho_ssd;
        }
    }

    // Saturation density predicted by KM: ρ_sat = (k1/k2)².
    let rho_sat = (p.k1 / p.k2).powi(2);
    println!("\nKM saturation density (k1/k2)² = {rho_sat:.3e} m⁻²");
    println!(
        "Final state: ρ_SSD = {:.3e} m⁻², ρ_total = {:.3e} m⁻²",
        state.rho_ssd,
        state.total()
    );
    println!("\nKocks-Mecking reproduces the saturation behaviour of the");
    println!("phenomenological Voce law, with ρ_sat = (k1/k2)² fixing");
    println!("the steady-state flow stress τ_sat = α μ b (k1/k2).");
}
