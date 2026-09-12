//! Golden regression test against `test-data/golden/taylor-factor/`.
//!
//! Recomputes the exact Bishop–Hill Taylor factors (classical axes + a
//! seeded random sweep) and compares against the committed baseline.
//! Regenerate deliberately with `cargo run -p golden-generate`.

use serde_json::Value;
use tpt_mat_crystallography::{CrystalStructure, SlipSystem};
use tpt_mat_rve::bishop_hill_taylor_factor_axis;

fn golden(path: &str) -> Value {
    let full = format!("{}/../../{}", env!("CARGO_MANIFEST_DIR"), path);
    let text = std::fs::read_to_string(&full)
        .unwrap_or_else(|e| panic!("missing golden file {full}: {e}"));
    serde_json::from_str(&text).unwrap()
}

fn close(name: &str, computed: f64, stored: f64) {
    let scale = stored.abs().max(1.0);
    assert!(
        (computed - stored).abs() <= 1.0e-9 * scale + 1.0e-12,
        "{name}: computed {computed:.16e} != golden {stored:.16e}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    *seed as f64 / u64::MAX as f64
}

/// Uniform random unit vector (Marsaglia sphere sampling) — must match
/// the generator exactly.
fn random_unit_axis(seed: &mut u64) -> [f64; 3] {
    loop {
        let u = 2.0 * lcg(seed) - 1.0;
        let v = 2.0 * lcg(seed) - 1.0;
        let s = u * u + v * v;
        if s < 1.0 && s > 0.0 {
            let f = 2.0 * (1.0 - s).sqrt();
            return [u * f, v * f, 1.0 - 2.0 * s];
        }
    }
}

fn check_file(path: &str, label: &str, slips: Option<&[SlipSystem]>) {
    let v = golden(path);
    for entry in v["axes"].as_array().unwrap() {
        let name = entry["name"].as_str().unwrap();
        let axis = [
            entry["axis"][0].as_f64().unwrap(),
            entry["axis"][1].as_f64().unwrap(),
            entry["axis"][2].as_f64().unwrap(),
        ];
        let m = bishop_hill_taylor_factor_axis(axis, 1.0, slips).taylor_factor;
        close(&format!("{label}:{name}"), m, entry["m"].as_f64().unwrap());
    }

    let n_total = v["n_random_total"].as_u64().unwrap() as usize;
    let mut seed = 0x5EED_1234_ABCD_0001u64;
    let mut sum = 0.0_f64;
    for i in 0..n_total {
        let axis = random_unit_axis(&mut seed);
        let m = bishop_hill_taylor_factor_axis(axis, 1.0, slips).taylor_factor;
        sum += m;
        if i < 64 {
            close(
                &format!("{label}:random[{i}]"),
                m,
                v["random_axes"][i]["m"].as_f64().unwrap(),
            );
        }
    }
    close(
        &format!("{label}:random_average"),
        sum / n_total as f64,
        v["random_average"].as_f64().unwrap(),
    );
    assert_eq!(
        n_total, 256,
        "golden n_random_total must stay in sync with the generator"
    );
}

#[test]
fn fcc_taylor_factor_golden() {
    check_file("test-data/golden/taylor-factor/fcc.json", "fcc", None);
}

#[test]
fn bcc_taylor_factor_golden() {
    let bcc = CrystalStructure::BCC.slip_systems();
    check_file("test-data/golden/taylor-factor/bcc.json", "bcc", Some(&bcc));
}
