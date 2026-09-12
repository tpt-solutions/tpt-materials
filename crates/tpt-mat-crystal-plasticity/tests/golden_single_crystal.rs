//! Golden regression test against `test-data/golden/
//! crystal-plasticity/fcc-single-crystal-tension.json`.
//!
//! Recomputes the 8-step uniaxial single-point sweep (Voce hardening,
//! FCC {111}<110>) and compares per-step stress / slip / hardening.

use serde_json::Value;

use tpt_mat_crystal_plasticity::{
    solve_increment_single_point, CrystalPlasticityModel, SymmetricFourthOrder,
};
use tpt_mat_crystallography::CrystalStructure;
use tpt_mat_hardening::{Hardening, HardeningState, VoceHardening, VoceParams};

fn golden() -> Value {
    let full = format!(
        "{}/../../test-data/golden/crystal-plasticity/fcc-single-crystal-tension.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&full)
        .unwrap_or_else(|e| panic!("missing golden file {full}: {e}"));
    serde_json::from_str(&text).unwrap()
}

fn close(name: &str, computed: f64, stored: f64) {
    let scale = stored.abs().max(1.0);
    assert!(
        (computed - stored).abs() <= 1.0e-9 * scale + 1.0e-6,
        "{name}: computed {computed:.8e} != golden {stored:.8e}"
    );
}

#[test]
fn fcc_single_crystal_tension_golden() {
    let v = golden();
    let elastic = SymmetricFourthOrder::cubic(168.4e9, 121.4e9, 75.4e9);
    let model = CrystalPlasticityModel::from_crystal_structure(
        CrystalStructure::FCC,
        Hardening::Voce(VoceHardening::uniform(VoceParams {
            tau_0: 30.0e6,
            tau_s: 60.0e6,
            theta_0: 500.0e6,
            gamma_c: 0.05,
        })),
        Default::default(),
        elastic,
    )
    .unwrap();
    let mut state = HardeningState::from_crss(&model.slip_systems);
    let steps = v["steps"].as_array().unwrap();
    for (i, step) in steps.iter().enumerate() {
        let k = i + 1;
        let eps_inc = tpt_math_linalg_fixed::Vec6::new(0.001 * k as f64, 0.0, 0.0, 0.0, 0.0, 0.0);
        let update = solve_increment_single_point(&model, &mut state, &eps_inc).unwrap();
        close(
            &format!("step[{k}].sigma_11"),
            update.stress.data[0],
            step["sigma_11"].as_f64().unwrap(),
        );
        close(
            &format!("step[{k}].sum_abs_slip_rates"),
            update.slip_rates.iter().map(|g| g.abs()).sum(),
            step["sum_abs_slip_rates"].as_f64().unwrap(),
        );
        close(
            &format!("step[{k}].sum_accumulated_shear"),
            state.accumulated_shear.iter().sum(),
            step["sum_accumulated_shear"].as_f64().unwrap(),
        );
        close(
            &format!("step[{k}].mean_crss"),
            state.crss.iter().sum::<f64>() / state.crss.len() as f64,
            step["mean_crss"].as_f64().unwrap(),
        );
    }
    // Physical sanity: monotone stress with plastic strain accumulation.
    let sigma = steps
        .iter()
        .map(|s| s["sigma_11"].as_f64().unwrap())
        .collect::<Vec<_>>();
    let shear = steps
        .iter()
        .map(|s| s["sum_accumulated_shear"].as_f64().unwrap())
        .collect::<Vec<_>>();
    assert!(sigma.windows(2).all(|w| w[1] > w[0]), "stress must grow");
    assert!(
        shear.windows(2).all(|w| w[1] >= w[0]),
        "accumulated shear must be monotone"
    );
}