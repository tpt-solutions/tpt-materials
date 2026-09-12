//! Golden regression test against `test-data/golden/
//! phase-field/spinodal-decomposition.json`.

use serde_json::Value;

use tpt_mat_phase_field::{BulkEnergy, PhaseFieldSolver, RegularSolutionParams};
use tpt_science::Grid2D;

fn golden() -> Value {
    let full = format!(
        "{}/../../test-data/golden/phase-field/spinodal-decomposition.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&full)
        .unwrap_or_else(|e| panic!("missing golden file {full}: {e}"));
    serde_json::from_str(&text).unwrap()
}

fn close(name: &str, computed: f64, stored: f64) {
    let scale = stored.abs().max(1.0);
    assert!(
        (computed - stored).abs() <= 1.0e-9 * scale + 1.0e-12,
        "{name}: computed {computed:.12e} != golden {stored:.12e}"
    );
}

#[test]
fn spinodal_decomposition_golden() {
    let v = golden();
    let n = v["grid_n"].as_u64().unwrap() as usize;
    let n_steps = v["n_steps"].as_u64().unwrap() as usize;
    let grid = Grid2D::new(n, n, 1.0);
    let mut c = vec![0.5_f64; n * n];
    for (idx, val) in c.iter_mut().enumerate() {
        let i = idx / n;
        let j = idx % n;
        let phi = 2.0 * std::f64::consts::PI * ((i as f64) * 0.05 + (j as f64) * 0.07);
        *val = 0.5 + 0.05 * phi.sin() + 0.03 * (2.0 * phi).cos();
    }
    let mut solver = PhaseFieldSolver::new_cahn_hilliard(
        grid,
        0.001,
        1.0,
        1.0,
        BulkEnergy::RegularSolution(RegularSolutionParams {
            omega: 4.0,
            rt: 1.0,
        }),
        &c,
    )
    .unwrap();
    let snap0 = solver.snapshot();
    close(
        "mean_initial",
        c.iter().sum::<f64>() / c.len() as f64,
        v["mean_initial"].as_f64().unwrap(),
    );
    solver.step_many(n_steps).unwrap();
    let snap = solver.snapshot();
    close(
        "energy_initial",
        snap0.free_energy,
        v["energy_initial"].as_f64().unwrap(),
    );
    close(
        "energy_final",
        snap.free_energy,
        v["energy_final"].as_f64().unwrap(),
    );
    close(
        "interface_area_initial",
        snap0.interface_area,
        v["interface_area_initial"].as_f64().unwrap(),
    );
    close(
        "interface_area_final",
        snap.interface_area,
        v["interface_area_final"].as_f64().unwrap(),
    );
    close(
        "mean_final",
        solver.concentration.iter().sum::<f64>() / solver.concentration.len() as f64,
        v["mean_final"].as_f64().unwrap(),
    );
    assert!(
        snap.free_energy <= snap0.free_energy + 1.0e-6,
        "free energy must not increase"
    );
}
