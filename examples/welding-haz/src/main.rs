//! Example: GMAW weld on AISI 4140 steel.
//!
//! 1. Heat-input / HAZ width from the Rosenthal line-source solution.
//! 2. CG-HAZ grain-coarsening + phase-fraction prediction.

use tpt_mat_welding::{
    heat_affected_zone, peak_temperature_profile, predict_haz_microstructure, BaseMetal,
    RosenthalWeld, WeldModel, WeldProcess,
};

fn main() {
    let weld = WeldModel::new(
        BaseMetal::aisi_4140(),
        None,
        WeldProcess::Gmaw,
        1500.0,
        0.005,
        0.8,
    );

    let haz = heat_affected_zone(&weld, weld.base_metal.a1_temperature);
    println!(
        "GMAW on AISI 4140 (Q = {:.0} J/m, v = {:.0} mm/s, η = {:.1})",
        weld.linear_heat_input,
        weld.scan_speed * 1000.0,
        weld.efficiency
    );
    println!("=============================================");
    println!(
        "HAZ width (above A₁ = {:.0} K): {:.2} mm",
        weld.base_metal.a1_temperature,
        haz.width * 1000.0
    );
    println!(
        "Peak temperature at fusion line: {:.0} K",
        haz.peak_temperature
    );

    let rw = RosenthalWeld {
        linear_heat_input: weld.linear_heat_input,
        scan_speed: weld.scan_speed,
        thermal_conductivity: weld.base_metal.thermal_conductivity,
        thermal_diffusivity: weld.base_metal.thermal_diffusivity,
        preheat_temperature: weld.base_metal.preheat_temperature,
    };
    let profile = peak_temperature_profile(&rw, 12);
    println!();
    println!("Peak-temperature profile:");
    for (y, t) in profile
        .positions
        .iter()
        .zip(profile.peak_temperatures.iter())
    {
        println!("  y = {y:.4e} m    T_peak = {t:>7.1} K");
    }

    let micro = predict_haz_microstructure(&weld, 100.0, 873.0);
    println!();
    println!("CG-HAZ prediction (100 s hold at 873 K):");
    println!("  CG-HAZ cooling rate = {:.1e} K/s", micro.cooling_rate);
    println!(
        "  CG-HAZ grain size   = {:.2} µm",
        micro.cg_haz_grain_size * 1.0e6
    );
    let labels = [
        "base_metal",
        "sub_critical",
        "inter_critical",
        "fine_grained",
        "coarse_grained",
        "martensite",
        "bainite",
        "pearlite",
    ];
    for (i, f) in micro.phase_fractions.iter().enumerate() {
        println!("  {label:<14} : {frac:.3}", label = labels[i], frac = f);
    }
}
