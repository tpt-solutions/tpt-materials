//! Heat-affected-zone (HAZ) from the Rosenthal line-source solution.
//!
//! The 2-D Rosenthal temperature field gives
//!
//! ```text
//! T(y) − T_0 = Q / (2 π k y) · exp(−v (y) / (2 α))
//! ```
//!
//! at the position of peak temperature `y` below the arc.  We
//! report the HAZ width `W_HAZ` as the perpendicular distance over
//! which the peak temperature exceeds `T_A1` (or some caller-
//! supplied threshold).

use serde::{Deserialize, Serialize};

use super::{WeldModel, WeldProcess};

/// Rosenthal welding solution inputs.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RosenthalWeld {
    /// Linear heat input `Q` (J/m).
    pub linear_heat_input: f64,
    /// Scan speed `v` (m/s).
    pub scan_speed: f64,
    /// Thermal conductivity `k` (W/m·K).
    pub thermal_conductivity: f64,
    /// Thermal diffusivity `α` (m²/s).
    pub thermal_diffusivity: f64,
    /// Pre-heat `T_0` (K).
    pub preheat_temperature: f64,
}

/// HAZ summary result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HazResult {
    /// HAZ width (m) — perpendicular distance over which `T_peak >
    /// T_threshold`.
    pub width: f64,
    /// Peak temperature at the fusion-line boundary (K).
    pub peak_temperature: f64,
    /// Threshold temperature used (K).
    pub threshold_temperature: f64,
}

/// Peak-temperature profile along a perpendicular to the weld line.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HazPeakProfile {
    /// Perpendicular positions (m).
    pub positions: Vec<f64>,
    /// Peak temperatures at each position (K).
    pub peak_temperatures: Vec<f64>,
}

/// Compute the HAZ width and peak temperature.
pub fn heat_affected_zone(weld: &WeldModel, threshold_temperature: f64) -> HazResult {
    let rw = RosenthalWeld {
        linear_heat_input: weld.linear_heat_input,
        scan_speed: weld.scan_speed,
        thermal_conductivity: weld.base_metal.thermal_conductivity,
        thermal_diffusivity: weld.base_metal.thermal_diffusivity,
        preheat_temperature: weld.base_metal.preheat_temperature,
    };
    // Bisect for the perpendicular distance y at which T_peak = T_thr.
    let mut lo = 1.0e-6;
    let mut hi = 1.0e-1;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let t = peak_temperature_at(&rw, mid);
        if t > threshold_temperature {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let width = lo;
    let peak_temperature = peak_temperature_at(&rw, width * 0.5);
    HazResult {
        width,
        peak_temperature,
        threshold_temperature,
    }
}

/// Peak temperature at perpendicular distance `y`.
pub fn peak_temperature_at(rw: &RosenthalWeld, y: f64) -> f64 {
    if y <= 0.0 {
        return f64::INFINITY;
    }
    let arg = -(rw.scan_speed * y) / (2.0 * rw.thermal_diffusivity);
    let exp_term = if arg < -50.0 { 0.0 } else { arg.exp() };
    rw.preheat_temperature
        + rw.linear_heat_input / (2.0 * core::f64::consts::PI * rw.thermal_conductivity * y)
            * exp_term
}

/// Peak-temperature profile sampled on a log-spaced perpendicular.
pub fn peak_temperature_profile(rw: &RosenthalWeld, n: usize) -> HazPeakProfile {
    let mut positions = Vec::with_capacity(n);
    let mut peak_temperatures = Vec::with_capacity(n);
    for i in 0..n {
        let frac = i as f64 / (n as f64 - 1.0).max(1.0);
        let log_y = -6.0 + frac * 4.0;
        let y = 10_f64.powf(log_y);
        positions.push(y);
        peak_temperatures.push(peak_temperature_at(rw, y));
    }
    HazPeakProfile {
        positions,
        peak_temperatures,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BaseMetal;

    fn weld() -> WeldModel {
        WeldModel::new(
            BaseMetal::aisi_4140(),
            None,
            WeldProcess::Gmaw,
            1000.0 * 25.0 / 0.005, // V=25V, I=200A, v=5 mm/s
            0.005,
            0.8,
        )
    }

    #[test]
    fn peak_temperature_decreases_with_distance() {
        let rw = RosenthalWeld {
            linear_heat_input: 2000.0 / 0.005,
            scan_speed: 0.005,
            thermal_conductivity: 42.0,
            thermal_diffusivity: 1.0e-5,
            preheat_temperature: 298.0,
        };
        let t1 = peak_temperature_at(&rw, 1.0e-3);
        let t2 = peak_temperature_at(&rw, 1.0e-2);
        assert!(t1 > t2);
    }

    #[test]
    fn haz_width_positive() {
        let h = heat_affected_zone(&weld(), 773.0);
        assert!(h.width > 0.0);
        assert!(h.width < 1.0);
    }

    #[test]
    fn profile_contains_n_samples() {
        let rw = RosenthalWeld {
            linear_heat_input: 2000.0 / 0.005,
            scan_speed: 0.005,
            thermal_conductivity: 42.0,
            thermal_diffusivity: 1.0e-5,
            preheat_temperature: 298.0,
        };
        let p = peak_temperature_profile(&rw, 50);
        assert_eq!(p.positions.len(), 50);
        assert_eq!(p.peak_temperatures.len(), 50);
    }
}
