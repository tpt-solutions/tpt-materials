//! Example: Al-Cu age-hardening curve.
//!
//! Demonstrates the precipitation-strengthening workflow by
//! sweeping mean precipitate radius (the time-evolution of the
//! population during ageing) and computing the Orowan bypass
//! strengthening increment at each stage.  Real age-hardening
//! curves couple KWN size-class evolution to LSW coarsening; the
//! functional form `Δτ(r̄, f)` is the focus here.
//!
//! Run with `cargo run -p precipitation-age-hardening`.

use tpt_mat_precipitation::{
    classical_nucleation_rate, lsw_coarsening_rate, nucleation_barrier,
    orowan_bypass_strengthening, shearing_strengthening,
};

const SHEAR_MODULUS_AL: f64 = 26.0e9;
const BURGERS_VECTOR_AL: f64 = 2.86e-10;
const POISSON_AL: f64 = 0.345;

fn main() {
    println!("Al-Cu age-hardening curve (Δτ vs r̄)");
    println!("=====================================");

    // Classical nucleation barrier and rate at 190 °C (463 K).
    let gamma = 0.05; // coherent nucleus (J/m²)
    let delta_gv = 5.0e9; // J/m³
    let delta_g_star = nucleation_barrier(gamma, delta_gv);
    let i_nuc = classical_nucleation_rate(1.0e40, 1.0e-20, 463.0, gamma, delta_gv);
    println!("Nucleation barrier ΔG* = {delta_g_star:.3e} J");
    println!("Initial nucleation rate I = {i_nuc:.3e} m⁻³·s⁻¹");

    // LSW coarsening rate at peak-aged conditions.
    let k_lsw = lsw_coarsening_rate(gamma, 1.0e-15, 0.005, 1.0e-5);
    println!("LSW coarsening rate constant K = {k_lsw:.3e} m³/s\n");

    // Typical ageing: Vf grows rapidly to ~0.04 during nucleation,
    // then is conserved as precipitates coarsen.
    let vf_peak = 0.04_f64;
    let vf_initial = 0.005_f64;

    println!("Ageing stages (Vf rises from {vf_initial} → {vf_peak}, then constant):\n");
    println!("  stage          r̄ (nm)    Vf        Δτ_Orowan (MPa)   Δτ_shear (MPa)");
    println!("  ------------   -------    -------   ---------------   --------------");
    let stages = [
        ("as-quenched", 0.0, 0.0),
        ("under-aged", 1.5e-9, vf_initial * 0.5),
        ("early-peak", 3.0e-9, vf_initial),
        ("peak-aged", 5.0e-9, vf_peak),
        ("peak-aged (θ')", 8.0e-9, vf_peak),
        ("over-aged", 15.0e-9, vf_peak),
        ("over-aged (coarse)", 30.0e-9, vf_peak),
        ("severely over-aged", 80.0e-9, vf_peak),
    ];
    for (label, r_bar, f_v) in &stages {
        let dt_orowan = orowan_bypass_strengthening(
            SHEAR_MODULUS_AL,
            BURGERS_VECTOR_AL,
            POISSON_AL,
            *r_bar,
            *f_v,
        );
        let dt_shear = shearing_strengthening(SHEAR_MODULUS_AL, BURGERS_VECTOR_AL, *r_bar, *f_v);
        println!(
            "  {:<14} {:>7.2}    {:>5.3}    {:>13.3}    {:>12.3}",
            label,
            r_bar * 1.0e9,
            f_v,
            dt_orowan / 1.0e6,
            dt_shear / 1.0e6,
        );
    }

    println!("\nInterpretation:");
    println!("  - At small r̄, shearing dominates (coherent particles).");
    println!("  - At large r̄, Orowan bypass dominates (incoherent particles).");
    println!("  - Peak hardness occurs near the shearing → Orowan crossover,");
    println!("    typically at r̄ ≈ 5–10 nm for Al-Cu θ'.");
}
