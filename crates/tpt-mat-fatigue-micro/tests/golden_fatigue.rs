//! Golden regression test against `test-data/golden/
//! degradation/fatigue-crack-initiation.json`.

use serde_json::Value;

use tpt_mat_fatigue_micro::{predict_crack_initiation, FatigueCriterion};

fn golden() -> Value {
    let full = format!(
        "{}/../../test-data/golden/degradation/fatigue-crack-initiation.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&full)
        .unwrap_or_else(|e| panic!("missing golden file {full}: {e}"));
    serde_json::from_str(&text).unwrap()
}

fn close(name: &str, computed: f64, stored: f64) {
    let scale = stored.abs().max(1.0);
    assert!(
        (computed - stored).abs() <= 1.0e-9 * scale + 1.0e-9,
        "{name}: computed {computed:.12e} != golden {stored:.12e}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    *seed as f64 / u64::MAX as f64
}

#[test]
fn fatigue_crack_initiation_golden() {
    let v = golden();
    let cases = v["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 3, "golden case count changed");

    // Rebuild the deterministic shear history used by the generator.
    let mut seed = 0xC0FFEE_2026u64;
    let mut shear_history: Vec<Vec<f64>> = (0..6)
        .map(|_| (0..12).map(|_| 0.001 + 0.009 * lcg(&mut seed)).collect())
        .collect();
    shear_history[3] = vec![0.04; 12];
    shear_history[4] = vec![0.0; 12];

    let criteria: [(&str, FatigueCriterion); 3] = [
        ("findley k=0.3", FatigueCriterion::Findley { k: 0.3 }),
        (
            "fatemi-socie",
            FatigueCriterion::FatemiSocie {
                sigma_y: 200.0,
                k: 0.8,
            },
        ),
        (
            "crystallographic gamma_c=0.5",
            FatigueCriterion::CrystallographicSlip {
                critical_accumulated_shear: 0.5,
            },
        ),
    ];

    for (idx, entry) in cases.iter().enumerate() {
        let expected_shear = entry["accumulated_shear"]
            .as_array()
            .unwrap()
            .iter()
            .map(|g| {
                g.as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_f64().unwrap())
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<_>>();
        assert_eq!(
            expected_shear.len(),
            shear_history.len(),
            "case[{idx}] grain count"
        );
        for (gi, (computed_g, stored_g)) in shear_history.iter().zip(&expected_shear).enumerate() {
            assert_eq!(
                computed_g.len(),
                stored_g.len(),
                "case[{idx}] grain {gi} system count"
            );
            for (si, (c, s)) in computed_g.iter().zip(stored_g).enumerate() {
                close(&format!("case[{idx}] history[{gi}][{si}]"), *c, *s);
            }
        }
        let (name, criterion) = &criteria[idx];
        let r = predict_crack_initiation(&shear_history, criterion.clone(), 1000.0, 0.25);
        close(
            &format!("{name}:cycles_to_initiation"),
            r.cycles_to_initiation,
            entry["cycles_to_initiation"].as_f64().unwrap(),
        );
        assert_eq!(
            r.critical_grain,
            entry["critical_grain"].as_u64().unwrap() as usize,
            "{name}: critical grain mismatch"
        );
        close(
            &format!("{name}:critical_fip"),
            r.critical_fip,
            entry["critical_fip"].as_f64().unwrap(),
        );
    }

    // Sanity: the highly sheared grain 3 dominates the FIP field.
    assert_eq!(
        cases[0]["critical_grain"].as_u64().unwrap(),
        3,
        "critical grain should be the heavily sheared grain"
    );
}
