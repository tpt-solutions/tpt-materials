//! Demo of the spec §11 backlog modules:
//!
//! 1. High-temperature oxidation (`tpt-mat-corrosion::oxidation`).
//! 2. Grain-boundary segregation + character distribution
//!    (`tpt-mat-grain-growth::interfaces`).
//! 3. JMAK recrystallization kinetics (`tpt-mat-grain-growth::recrystallization`).
//! 4. Stereological grain-size / phase-fraction quantification
//!    (`tpt-mat-grain-growth::stereology`).
//! 5. Closed-form inverse calibration (`tpt-mat-inverse`).
//!
//! See RFC 0016 (`rfcs/0016-backlog-modules.md`) for the underlying
//! modelling assumptions and references.

use tpt_mat_corrosion::{
    BreakawayCriterion, LinearBreakawayRate, OxidationModel, ParabolicRateConstant, ScaleGrowth,
};
use tpt_mat_grain_growth::{
    AreaFraction2D, GBCDEntry, GrainBoundaryCharacterDistribution, GrainBoundaryEnergy,
    JmakRecrystallization, LangmuirMcLean, LinearIntercept, TripleJunction, ZenerHollomon,
};
use tpt_mat_inverse::{arrhenius_fit, norton_creep_fit, power_law_sn_fit};

const R: f64 = 8.314_462_618;

fn oxidation_demo() {
    println!("=== High-temperature oxidation ===");
    let kp = ParabolicRateConstant {
        k_p_0: 1.0e-6,
        activation_energy: 150.0e3,
    };
    let model = OxidationModel {
        parabolic: ScaleGrowth::at_temperature(kp, 900.0, R),
        breakaway: BreakawayCriterion {
            critical_thickness: 5.0e-6,
            critical_strain: 1.0,
        },
        linear: LinearBreakawayRate::new(1.0e-9),
    };
    for &(t, label) in &[
        (60.0_f64, "1 min"),
        (3600.0, "1 h"),
        (86400.0, "1 d"),
        (86400.0 * 365.25, "1 yr"),
    ] {
        let x = model.thickness(t, 0.0);
        println!("  t = {:>10}: scale thickness = {:.2e} m", label, x);
    }
}

fn interfaces_demo() {
    println!("\n=== Grain-boundary segregation (Langmuir–McLean) ===");
    let seg = LangmuirMcLean::new(-35.0e3);
    for &t in &[600.0_f64, 800.0, 1000.0] {
        let beta = seg.enrichment_factor(0.01, t, R);
        println!(
            "  T = {:>4.0} K, X_bulk = 0.01 -> X_GB enrichment β = {:.2}",
            t, beta
        );
    }

    println!("\n=== Triple-junction force balance ===");
    let tj = TripleJunction {
        gammas: [1.0, 1.2, 0.9],
    };
    for i in 0..3 {
        if let Some(phi) = tj.dihedral_angle(i) {
            println!("  φ_{} = {:.2}°", i, phi.to_degrees());
        }
    }

    println!("\n=== GB character distribution ===");
    let mut gbcd = GrainBoundaryCharacterDistribution::new();
    for (theta_deg, weight) in [
        (3.0_f64, 5.0_f64),
        (8.0, 8.0),
        (15.0, 4.0),
        (60.0, 2.0),
        (38.0, 3.0),
    ] {
        gbcd.push(GBCDEntry {
            misorientation: theta_deg.to_radians(),
            plane_normal: [0.0, 0.0, 1.0],
            weight,
        });
    }
    let lab_frac = gbcd.low_angle_fraction(10.0_f64.to_radians());
    let csl_frac = gbcd.csl_fraction(9);
    println!(
        "  LAB fraction = {:.2}, special-boundary (Σ ≤ 9) fraction = {:.2}",
        lab_frac, csl_frac
    );

    println!("\n=== Read–Shockley GB energy ===");
    let e = GrainBoundaryEnergy::default();
    for &theta_deg in &[2.0_f64, 5.0, 10.0, 30.0, 45.0] {
        let g = e.at_misorientation(theta_deg.to_radians());
        println!("  θ = {:>4.1}°  → γ_GB = {:.3} J/m²", theta_deg, g);
    }
}

fn recrystallization_demo() {
    println!("\n=== JMAK static recrystallization ===");
    let jmak = JmakRecrystallization::default();
    for &t in &[1.0_f64, 10.0, 100.0, 1000.0, 10000.0] {
        let x = jmak.recrystallized_fraction(t, 900.0, R);
        let t50 = jmak.t_50(900.0, R);
        println!(
            "  t = {:>6.0} s at 900 K: X = {:.3}, t_50 = {:.1} s",
            t, x, t50
        );
    }

    println!("\n=== Zener–Hollomon ===");
    let mut zh = ZenerHollomon::default();
    zh.strain_rate = 1.0;
    let sigma = zh.steady_state_stress(1000.0, R);
    println!(
        "  ε̇ = 1 s⁻¹, T = 1000 K: Z = {:.2e}, σ_s = {:.2} MPa",
        zh.z(1000.0, R),
        sigma
    );
}

fn stereology_demo() {
    println!("\n=== Stereology ===");
    let li = LinearIntercept::new(1.0e-3, 200);
    println!(
        "  Mean intercept = {:.2e} m, ASTM G = {:.1}",
        li.mean_intercept(),
        li.astm_grain_size_number()
    );
    let af = AreaFraction2D::new(312, 1024);
    println!(
        "  Phase area fraction = {:.3} ± {:.4}",
        af.area_fraction(),
        af.uncertainty()
    );
}

fn inverse_demo() {
    println!("\n=== Closed-form inverse calibration ===");

    // Synthesise an Arrhenius dataset.
    let d_0_true = 1.0e-5_f64;
    let q_true = 100.0e3_f64;
    let arr_data: Vec<(f64, f64)> = (800..1100).step_by(20).fold(Vec::new(), |mut v, t| {
        let t = t as f64;
        let ln_d = (d_0_true * (-q_true / (R * t)).exp()).ln();
        v.push((1.0 / t, ln_d));
        v
    });
    let fit = arrhenius_fit(&arr_data, R);
    println!(
        "  Arrhenius: D_0 recovered = {:.3e} (true {:.1e}), Q = {:.1} kJ/mol (true {:.0})",
        fit.k_0,
        d_0_true,
        fit.q / 1000.0,
        q_true / 1000.0
    );

    // Norton creep.
    let n_true = 4.0_f64;
    let a_true = 1.0e-5_f64;
    let sn: Vec<(f64, f64)> = (10..100).step_by(10).fold(Vec::new(), |mut v, s| {
        let s = s as f64;
        v.push((s.ln(), (a_true * s.powf(n_true)).ln()));
        v
    });
    let nf = norton_creep_fit(&sn);
    println!(
        "  Norton–Bailey: A = {:.3e}, n = {:.3} (true {:.0})",
        nf.a, nf.n, n_true
    );

    // Power-law S–N.
    let sigma_f_true = 1000.0_f64;
    let b_true = -0.08_f64;
    let sn_data: Vec<(f64, f64)> = (1..6).fold(Vec::new(), |mut v, i| {
        let n = 10.0_f64.powi(i);
        let s = sigma_f_true * (2.0 * n).powf(b_true);
        v.push((n, s));
        v
    });
    let pl = power_law_sn_fit(&sn_data);
    println!(
        "  Basquin S–N: σ_f' = {:.1} MPa (true {:.0}), b = {:.4} (true {:.2})",
        pl.sigma_f_prime, sigma_f_true, pl.b, b_true
    );
}

fn main() {
    oxidation_demo();
    interfaces_demo();
    recrystallization_demo();
    stereology_demo();
    inverse_demo();
}