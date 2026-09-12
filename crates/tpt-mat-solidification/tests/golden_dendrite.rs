//! Golden regression test against `test-data/golden/
//! phase-field/dendritic-solidification.json`.

use serde_json::Value;

use tpt_mat_phase_field::BulkEnergy;
use tpt_mat_solidification::{AnisotropyModel, AnisotropyMode, SolidificationSolver};
use tpt_science::Grid2D;

fn golden() -> Value {
    let full = format!(
        "{}/../../test-data/golden/phase-field/dendritic-solidification.json",
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
fn dendritic_solidification_golden() {
    let v = golden();
    let n = v["grid_n"].as_u64().unwrap() as usize;
    let n_steps = v["n_steps"].as_u64().unwrap() as usize;
    let grid = Grid2D::new(n, n, 1.0);
    let mut eta = vec![0.0_f64; n * n];
    let cx = n as f64 / 2.0;
    let cy = n as f64 / 2.0;
    let seed_radius = 4.0_f64;
    for i in 0..n {
        for j in 0..n {
            let d = ((i as f64 - cy).powi(2) + (j as f64 - cx).powi(2)).sqrt();
            eta[i * n + j] = if d < seed_radius { 1.0 } else { 0.0 };
        }
    }
    let anisotropy = AnisotropyModel {
        strength: 0.04,
        mode: AnisotropyMode::Cubic4Fold,
    };
    let mut solver = SolidificationSolver::new(
        grid,
        &eta,
        anisotropy,
        BulkEnergy::DoubleWell { well_depth: 1.0 },
        0.005,
        1.0,
        1.0,
    )
    .unwrap();
    let snap0 = solver.snapshot();
    for _ in 0..n_steps {
        solver.step().unwrap();
    }
    let snap = solver.snapshot();
    close(
        "solid_fraction_initial",
        snap0.solid_fraction,
        v["solid_fraction_initial"].as_f64().unwrap(),
    );
    close(
        "solid_fraction_final",
        snap.solid_fraction,
        v["solid_fraction_final"].as_f64().unwrap(),
    );
    close("tip_velocity", snap.tip_velocity, v["tip_velocity"].as_f64().unwrap());
    close(
        "primary_arm_spacing",
        snap.primary_arm_spacing,
        v["primary_arm_spacing"].as_f64().unwrap(),
    );
    close(
        "secondary_arm_spacing",
        snap.secondary_arm_spacing,
        v["secondary_arm_spacing"].as_f64().unwrap(),
    );
    assert!(
        snap.solid_fraction > snap0.solid_fraction,
        "solid fraction must grow"
    );
}