//! Golden regression test against `test-data/golden/
//! energy-materials/corrosion-polarization.json`.

use serde_json::Value;

use tpt_mat_corrosion::{
    corrosion_rate, polarization_curve, CorrosionModel, ElectrodeKinetics, PolarizationBranch,
};

fn golden() -> Value {
    let full = format!(
        "{}/../../test-data/golden/energy-materials/corrosion-polarization.json",
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

fn electrode_from(v: &Value) -> ElectrodeKinetics {
    ElectrodeKinetics::from_alphas(
        v["equilibrium_potential"].as_f64().unwrap(),
        v["exchange_current_density"].as_f64().unwrap(),
        v["alpha_a"].as_f64().unwrap(),
        v["alpha_c"].as_f64().unwrap(),
        v["n"].as_f64().unwrap(),
        298.0,
    )
}

#[test]
fn corrosion_polarization_golden() {
    let v = golden();
    let fe_model = CorrosionModel::new(
        electrode_from(&v["anode"]),
        electrode_from(&v["cathode"]),
        v["ph"].as_f64().unwrap(),
        v["temperature_k"].as_f64().unwrap(),
    );

    // Corrosion rate: Fe and Ti against the same H2 cathode.
    let cr_fe = corrosion_rate(&fe_model, 0.055_845, 2.0, 7874.0);
    let golden_fe = &v["corrosion_fe"];
    close(
        "fe:current_density",
        cr_fe.current_density,
        golden_fe["current_density"].as_f64().unwrap(),
    );
    close(
        "fe:corrosion_potential",
        cr_fe.corrosion_potential,
        golden_fe["corrosion_potential"].as_f64().unwrap(),
    );
    close(
        "fe:penetration_rate_mm_per_yr",
        cr_fe.penetration_rate_mm_per_yr,
        golden_fe["penetration_rate_mm_per_yr"].as_f64().unwrap(),
    );
    close(
        "fe:mass_loss_rate_g_per_m2_day",
        cr_fe.mass_loss_rate_g_per_m2_day,
        golden_fe["mass_loss_rate_g_per_m2_day"].as_f64().unwrap(),
    );

    let ti_anode = ElectrodeKinetics::from_alphas(-0.86, 1.0e-7, 0.5, 0.5, 3.0, 298.0);
    let ti_model = CorrosionModel::new(ti_anode, electrode_from(&v["cathode"]), 0.0, 298.0);
    let cr_ti = corrosion_rate(&ti_model, 0.047_867, 3.0, 4506.0);
    let golden_ti = &v["corrosion_ti"];
    close(
        "ti:current_density",
        cr_ti.current_density,
        golden_ti["current_density"].as_f64().unwrap(),
    );
    close(
        "ti:penetration_rate_mm_per_yr",
        cr_ti.penetration_rate_mm_per_yr,
        golden_ti["penetration_rate_mm_per_yr"].as_f64().unwrap(),
    );

    // Polarization curve (Net branch) over the golden potential range.
    let pol = &v["polarization"];
    let e_min = pol["e_min"].as_f64().unwrap();
    let e_max = pol["e_max"].as_f64().unwrap();
    let num_points = pol["num_points"].as_u64().unwrap() as usize;
    let pc = polarization_curve(
        &fe_model,
        (e_min, e_max),
        num_points,
        PolarizationBranch::Net,
    );
    let golden_potentials = pol["potentials"].as_array().unwrap();
    let golden_currents = pol["currents"].as_array().unwrap();
    assert_eq!(
        pc.potentials.len(),
        num_points,
        "polarization point count changed"
    );
    assert_eq!(golden_potentials.len(), num_points);
    assert_eq!(golden_currents.len(), num_points);
    for i in 0..num_points {
        close(
            &format!("polarization:potential[{i}]"),
            pc.potentials[i],
            golden_potentials[i].as_f64().unwrap(),
        );
        close(
            &format!("polarization:current[{i}]"),
            pc.currents[i],
            golden_currents[i].as_f64().unwrap(),
        );
    }

    // Sanity: the Fe Net polarization curve changes sign across E_corr,
    // i.e. the computed mixed potential is where anode == cathode.
    let e_corr = golden_fe["corrosion_potential"].as_f64().unwrap();
    let crossed = pc
        .potentials
        .windows(2)
        .zip(pc.currents.windows(2))
        .any(|(e, i)| e[0] <= e_corr && e_corr < e[1] && i[0] * i[1] <= 0.0);
    assert!(
        crossed,
        "Fe polarization curve must cross zero near E_corr = {e_corr}"
    );
}
