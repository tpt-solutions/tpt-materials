//! Example: hydrogen accumulation at a notch-tip stress field.
//!
//! Models the steady-state hydrogen concentration around a
//! blunt notch under remote tension.  The hydrostatic stress
//! gradient drives an uphill flux via
//! `J = -D C V_H / (RT) ∇σ_h`.  At equilibrium the lattice
//! concentration follows `C ∝ exp(V_H σ_h / RT)`.
//!
//! Outputs:
//! 1. Hydrostatic stress profile `σ_h(x)` ahead of the notch.
//! 2. Hydrogen concentration enrichment profile `C(x)/C_∞`.
//! 3. HEDE and HELP susceptibility at the notch root.
//!
//! Run with `cargo run -p hydrogen-embrittlement-notch`.

use tpt_mat_hydrogen_embrittlement::{
    hede_threshold, help_threshold, hydrogen_diffusivity_effective, stress_driven_flux,
    susceptibility_index, HedeParams, HelpParams, OrianiParams,
};

const R_GAS: f64 = 8.314_462;
const V_H_M3_PER_MOL: f64 = 2.0e-6; // H partial-molar volume in steel

fn hydrostatic_stress_ahead_of_notch(x: f64, notch_radius: f64, q_max: f64) -> f64 {
    // Simple elastic approximation: σ_h peaks at the notch root
    // and decays ∝ 1 / (1 + x/ρ)².  q_max is the peak value.
    let decay = 1.0 / (1.0 + x / notch_radius).powi(2);
    q_max * decay
}

/// Triaxiality ratio ahead of the notch (σ_m / σ_eff).
fn triaxiality(x: f64, notch_radius: f64) -> f64 {
    // For a notch in elastic / elastic-plastic regime,
    // triaxiality starts at ~0.5 at the root and approaches
    // 1/3 far from it.
    let ramp = (1.0 - (-x / notch_radius).exp()).clamp(0.0, 1.0);
    0.5 + (1.0 / 3.0 - 0.5) * ramp
}

fn main() {
    let c_infinity = 1.0e-3_f64; // mol/m³ bulk H concentration
    let q_max = 800.0e6_f64; // peak hydrostatic stress (Pa)
    let notch_radius = 50.0e-6_f64; // 50 µm
    let t_k = 300.0_f64;

    let oriani = OrianiParams {
        d_lattice: 1.0e-8,
        trap_density: 1.0e25,
        binding_energy_j_per_mol: 30_000.0,
        temperature_k: t_k,
    };
    let d_eff = hydrogen_diffusivity_effective(c_infinity, &oriani);
    let d_lattice = oriani.d_lattice;
    println!("Hydrogen transport at a notch-tip (V_H = {V_H_M3_PER_MOL:.1e} m³/mol, T = {t_k} K)");
    println!("====================================================================");
    println!("  D_L     = {d_lattice:.3e} m²/s");
    println!("  D_eff   = {d_eff:.3e} m²/s (Oriani local equilibrium)");
    println!("  ratio   = {:.3}", d_eff / d_lattice);

    println!("\n  x (µm)    σ_h (MPa)    C(x)/C∞    C(x) (mol/m³)    Susc.    HEDE    HELP");
    println!("  -------   ---------    -------    -------------    ------   -----   ----");

    let x_values_um = [0.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0, 500.0, 1000.0];
    let hede = HedeParams { c_critical: 5.0e-3 };
    let help = HelpParams {
        c_critical: 2.0e-3,
        triaxiality_critical: 0.4,
    };

    for &x_um in &x_values_um {
        let x_m = x_um * 1.0e-6;
        let sigma_h = hydrostatic_stress_ahead_of_notch(x_m, notch_radius, q_max);
        // Equilibrium concentration under stress (Sievert's-law-like):
        let enrichment = (V_H_M3_PER_MOL * sigma_h / (R_GAS * t_k)).exp();
        let c_local = c_infinity * enrichment;
        let triax = triaxiality(x_m, notch_radius);
        let susc = susceptibility_index(c_local, triax, &hede, &help);
        let hede_on = hede_threshold(c_local, &hede);
        let help_on = help_threshold(c_local, triax, &help);
        // Show local flux as sanity (∂σ_h/∂x).
        let d_sigma_h = -2.0 * q_max / (notch_radius * (1.0 + x_m / notch_radius).powi(3));
        let _flux = stress_driven_flux(d_eff, c_local, V_H_M3_PER_MOL, t_k, d_sigma_h);

        println!(
            "  {x_um:>6.0}    {sigma_h:>9.3}    {enrichment:>5.2}    {c_local:>11.3e}    {susc:>4.2}    {:>5}    {:>4}",
            if hede_on { "YES" } else { "no" },
            if help_on { "YES" } else { "no" },
        );
    }

    println!("\nInterpretation:");
    println!("  - Hydrostatic tension concentrates H at the notch tip (Sieverts enrichment).");
    println!("  - HEDE triggers once local C crosses the decohesion threshold.");
    println!("  - HELP additionally requires high triaxiality, often met at notch roots.");
}
