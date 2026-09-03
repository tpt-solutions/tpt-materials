//! Example: JMAK / Avrami TTT diagram + Koistinen–Marburger martensite.
//!
//! Builds a Time-Temperature-Transformation (TTT) diagram for a
//! hypothetical diffusional phase transformation, runs an isothermal hold
//! through the nose of the C curve, then performs a continuous-cooling
//! transformation (CCT) through the martensite-start temperature.
//!
//! Run with: `cargo run --release -p phase-transform-jmak`.

use tpt_mat_phase_transform::{
    AvramiModel, AvramiParams, KMParams, KoistinenMarburger, ThermalHistoryStep,
    TransformationSolver,
};

fn main() {
    // Diffusional transformation (e.g.g. austenite → pearlite).
    let avrami = AvramiModel::new(AvramiParams {
        k0: 1.0e-3,
        activation_energy_kj_per_mol: 60.0,
        r_kj_per_mol_k: 8.314_462_618e-3,
        avrami_exponent: 3.0,
    });

    // Diffusionless martensitic transformation with M_s = 600 K.
    let km = KoistinenMarburger::new(KMParams {
        m_start_k: 600.0,
        alpha_per_k: 0.011,
        retained_at_m_start: 0.01,
    });

    println!("=== TTT diagram (Avrami): time to 50 % transformed ===");
    println!("T (K)   t (s)");
    for tk in (600..=1100).step_by(50) {
        let t = avrami
            .time_to_fraction(0.5, tk as f64)
            .unwrap_or(f64::INFINITY);
        println!("{tk:>5}   {t:>10.3}");
    }

    println!("\n=== Isothermal hold at 900 K for 600 s ===");
    let mut solver = TransformationSolver::new(avrami, km, 1100.0);
    let checkpoints = [0.0, 100.0, 200.0, 400.0, 600.0];
    let mut last = 0.0;
    for &target in &checkpoints {
        let dt = target - last;
        if dt > 0.0 {
            solver.apply_isothermal(900.0, dt);
        }
        last = target;
        println!(
            "t = {:>5.0} s  f_av = {:.4}  f_km = {:.4}  f_total = {:.4}",
            solver.state.time_s,
            solver.state.avrami_fraction,
            solver.state.martensite_fraction,
            solver.state.combined_fraction
        );
    }

    println!("\n=== Continuous cooling from 1100 K to 200 K over 2000 s ===");
    let mut solver2 = TransformationSolver::new(avrami, km, 1100.0);
    solver2
        .apply_history(&[ThermalHistoryStep {
            t_initial_k: 1100.0,
            t_final_k: 200.0,
            duration_seconds: 2000.0,
        }])
        .unwrap();
    println!(
        "Final: f_av = {:.4}  f_km = {:.4}  f_total = {:.4}",
        solver2.state.avrami_fraction,
        solver2.state.martensite_fraction,
        solver2.state.combined_fraction
    );

    // Emit a simple JSON TTT curve.
    let mut json = String::from("{\"ttt\":[");
    for (i, tk) in (600..=1100).step_by(50).enumerate() {
        let t = avrami
            .time_to_fraction(0.5, tk as f64)
            .unwrap_or(f64::INFINITY);
        if i > 0 {
            json.push(',');
        }
        json.push_str(&format!("{{\"t_k\":{tk},\"t_s\":{t}}}"));
    }
    json.push_str("]}");
    eprintln!("\nTTT JSON: {}", json);
}
