//! Golden regression test against `test-data/golden/
//! energy-materials/battery-degradation.json`.

use serde_json::Value;

use tpt_mat_battery::{capacity_fade_curve, default_active_material, BatteryChemistry};

fn golden() -> Value {
    let full = format!(
        "{}/../../test-data/golden/energy-materials/battery-degradation.json",
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
fn battery_degradation_golden() {
    let v = golden();
    let material = default_active_material(BatteryChemistry::Nmc811);
    let num_cycles = v["num_cycles"].as_u64().unwrap() as usize;
    let curve = capacity_fade_curve(&material, 1.0, num_cycles, 298.0);
    close(
        "final_capacity_retention",
        *curve.capacity_retention.last().unwrap(),
        v["final_capacity_retention"].as_f64().unwrap(),
    );
    close(
        "final_resistance_growth",
        *curve.resistance_growth.last().unwrap(),
        v["final_resistance_growth"].as_f64().unwrap(),
    );
    for milestone in v["milestones"].as_array().unwrap() {
        let cycle = milestone["cycle"].as_u64().unwrap() as usize;
        close(
            &format!("milestone[{cycle}].capacity_retention"),
            curve.capacity_retention[cycle - 1],
            milestone["capacity_retention"].as_f64().unwrap(),
        );
        close(
            &format!("milestone[{cycle}].resistance_growth"),
            curve.resistance_growth[cycle - 1],
            milestone["resistance_growth"].as_f64().unwrap(),
        );
    }
    // Sanity: retention decreases and resistance grows with cycling.
    let cap = v["final_capacity_retention"].as_f64().unwrap();
    let res = v["final_resistance_growth"].as_f64().unwrap();
    assert!(cap > 0.0 && cap <= 1.0, "retention in (0, 1]");
    assert!(res > 0.0, "resistance must grow");
}
