//! Golden regression test against `test-data/golden/
//! degradation/gtn-void-growth.json`.

use serde_json::Value;

use tpt_mat_damage::{gtn_yield_function, update_porosity, GtnParams};

fn golden() -> Value {
    let full = format!(
        "{}/../../test-data/golden/degradation/gtn-void-growth.json",
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
fn gtn_void_growth_golden() {
    let v = golden();
    let params = GtnParams::classic();
    let mut f = v["f_initial"].as_f64().unwrap();
    for (i, stored) in v["f_trajectory"].as_array().unwrap().iter().enumerate() {
        f = update_porosity(f, 0.01, 0.05, &params);
        close(
            &format!("f_trajectory[{i}]"),
            f,
            stored.as_f64().unwrap(),
        );
    }
    for point in v["yield_grid"].as_array().unwrap() {
        let seq = point["sigma_eq_over_y"].as_f64().unwrap();
        let sh = point["sigma_h_over_y"].as_f64().unwrap();
        let fv = point["f"].as_f64().unwrap();
        let computed = gtn_yield_function(seq, sh, 1.0, fv, &params);
        close(
            &format!("yield(seq={seq}, sh={sh}, f={fv})"),
            computed,
            point["value"].as_f64().unwrap(),
        );
    }
    // Sanity: porosity grows monotonically and never exceeds the
    // coalescence/failure thresholds by an unbounded amount.
    let traj: Vec<f64> = v["f_trajectory"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    assert!(traj.windows(2).all(|w| w[1] > w[0]), "porosity must grow");
    assert!(traj.last().unwrap() < 1.0, "porosity must stay < 1");
}