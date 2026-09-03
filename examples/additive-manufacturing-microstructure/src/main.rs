//! Example: LPBF Ti-6Al-4V scan.
//!
//! 1. Compute the Rosenthal moving-point-source thermal history at
//!    a point 0.1 mm from the scan track.
//! 2. Run the Hunt columnar-to-equiaxed map at the cooling-rate point.
//! 3. Build a 1-D residual-stress field along the build direction.

use tpt_mat_additive::{
    cooling_rate, predict_microstructure, residual_stress_field, thermal_history, AmProcess,
    MovingPointSource, SolidificationMap, ThermalContraction,
};

fn main() {
    let process = AmProcess::LaserPowderBedFusion {
        laser_power: 200.0,
        scan_speed: 1.0, // 1 m/s
        hatch_spacing: 100.0e-6,
        layer_thickness: 30.0e-6,
        beam_radius: 50.0e-6,
    };
    println!("LPBF Ti-6Al-4V");
    println!("==============");
    println!(
        "  linear heat input     = {:.1} J/m",
        process.linear_heat_input()
    );
    println!(
        "  volumetric energy     = {:.3e} J/m³",
        process.volumetric_energy_density()
    );

    let ti_src = MovingPointSource {
        linear_heat_input: 200.0,
        scan_speed: 1.0,
        thermal_conductivity: 21.0, // Ti-6Al-4V
        thermal_diffusivity: 9.1e-6,
        preheat_temperature: 873.0, // build-plate pre-heat
    };
    let d = 100.0e-6;
    let cr = cooling_rate(&ti_src, d);
    println!();
    println!("Rosenthal at d = {d:.1e} m:");
    println!("  peak temperature   = {:.1} K", cr.peak_temperature);
    println!("  cooling rate       = {:.1e} K/s", cr.cooling_rate);

    let map = SolidificationMap::default();
    let micro = predict_microstructure(&cr, &map);
    println!();
    println!("Hunt CET map:");
    println!("  grain size         = {:.2} µm", micro.grain_size * 1.0e6);
    println!("  columnar fraction  = {:.3}", micro.phase_fractions[0]);
    println!("  equiaxed fraction  = {:.3}", micro.phase_fractions[1]);

    let tc = ThermalContraction::default();
    let rsf = residual_stress_field(&tc, 5.0e-3, 6);
    println!();
    println!("Residual-stress field:");
    for (x, s) in rsf.positions.iter().zip(rsf.stresses.iter()) {
        println!("  x = {x:.4e} m    σ = {s:.3e} Pa", x = x, s = s);
    }

    let h = thermal_history(&ti_src, d, 21);
    println!();
    println!("Sample thermal history ({} points):", h.times.len());
    for i in (0..h.times.len()).step_by(4) {
        println!(
            "  t = {:>7.4} s    T = {:>7.1} K",
            h.times[i], h.temperatures[i]
        );
    }
    println!(
        "  peak temperature reached: {:.1} K at t = {:.4} s",
        h.peak_temperature, h.time_at_peak
    );
}
