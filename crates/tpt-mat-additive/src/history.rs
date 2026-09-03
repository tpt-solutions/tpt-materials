//! Moving-source analytic thermal-history model.
//!
//! The 2-D Rosenthal (1941) solution for a moving point source on
//! a semi-infinite plate gives the temperature field
//!
//! ```text
//! T(r) − T_0 = Q / (2 π k) · exp(−v (ξ + R) / (2 α))
//!             / R
//! ```
//!
//! where `R = √(ξ² + y² + z²)` and `ξ = x − v t`.
//!
//! The 3-D Gaussian-beam correction (Eagar–Tsai 1983) uses an
//! effective radius `r_0` in place of a point source.

use serde::{Deserialize, Serialize};

/// A moving point source (Rosenthal).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MovingPointSource {
    /// Linear heat input `Q = P / v` (J/m).
    pub linear_heat_input: f64,
    /// Scan speed `v` (m/s).
    pub scan_speed: f64,
    /// Thermal conductivity `k` (W/m·K).
    pub thermal_conductivity: f64,
    /// Thermal diffusivity `α = k / (ρ c_p)` (m²/s).
    pub thermal_diffusivity: f64,
    /// Pre-heat / substrate temperature `T_0` (K).
    pub preheat_temperature: f64,
}

/// A Gaussian-beam heat source (Eagar–Tsai).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GaussianBeam {
    /// Total incident power `P` (W).
    pub power: f64,
    /// Beam 1/e² radius `r_0` (m).
    pub beam_radius: f64,
    /// Absorptivity `η` (0–1).
    pub absorptivity: f64,
    /// Scan speed `v` (m/s).
    pub scan_speed: f64,
    /// Thermal conductivity `k` (W/m·K).
    pub thermal_conductivity: f64,
    /// Thermal diffusivity `α` (m²/s).
    pub thermal_diffusivity: f64,
    /// Pre-heat temperature `T_0` (K).
    pub preheat_temperature: f64,
}

/// Cooling-rate sample (K/s).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CoolingRate {
    /// Peak temperature reached (K).
    pub peak_temperature: f64,
    /// Time at peak (s).
    pub time_at_peak: f64,
    /// Cooling rate at peak (K/s, positive on cooling).
    pub cooling_rate: f64,
    /// Distance from the scan track (m).
    pub distance: f64,
}

/// Sampled thermal history at a point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ThermalHistory {
    /// Time samples (s).
    pub times: Vec<f64>,
    /// Temperature samples (K).
    pub temperatures: Vec<f64>,
    /// Peak temperature (K).
    pub peak_temperature: f64,
    /// Time of peak (s).
    pub time_at_peak: f64,
    /// Estimated cooling rate at peak (K/s).
    pub cooling_rate_at_peak: f64,
}

/// Compute the moving-point-source temperature field at a sample
/// time.  `ξ = x − v t` (along-track), `r = √(ξ² + y²)`.
pub fn rosenthal_2d(src: &MovingPointSource, xi: f64, y: f64, z: f64) -> f64 {
    let r_sq = xi * xi + y * y + z * z;
    let r = r_sq.sqrt();
    if r < 1.0e-12 {
        return f64::INFINITY;
    }
    let v = src.scan_speed;
    let arg = -(v * (xi + r)) / (2.0 * src.thermal_diffusivity);
    let exp_term = if arg < -50.0 { 0.0 } else { arg.exp() };
    src.preheat_temperature
        + src.linear_heat_input / (2.0 * core::f64::consts::PI * src.thermal_conductivity * r)
            * exp_term
}

/// Compute the cooling rate and peak temperature at distance `d` from
/// the scan track using the Rosenthal solution.  Distance `d` is the
/// perpendicular distance (m).
pub fn cooling_rate(src: &MovingPointSource, distance: f64) -> CoolingRate {
    // The peak temperature at (0, d, 0) occurs when the source is
    // closest, i.e. when ξ = 0.
    let t_peak = 0.0;
    let t_peak_val = src.preheat_temperature
        + src.linear_heat_input
            / (2.0 * core::f64::consts::PI * src.thermal_conductivity * distance)
            * (-src.scan_speed * distance / (2.0 * src.thermal_diffusivity)).exp();
    // Empirical cooling-rate estimate from the line-source solution:
    // cooling ∝ v · G (scan speed × interface gradient).  We use the
    // cooling rate at the trailing edge of the heat-affected zone as a
    // proxy:
    //   cooling ≈ v · (T_peak − T_0) / (2 r_c)
    // where r_c is the characteristic thermal radius ~ 2α / v.
    let r_c = 2.0 * src.thermal_diffusivity / src.scan_speed.max(1.0e-9);
    let cooling = if r_c > 0.0 {
        src.scan_speed * (t_peak_val - src.preheat_temperature) / (2.0 * r_c)
    } else {
        0.0
    };
    CoolingRate {
        peak_temperature: t_peak_val,
        time_at_peak: t_peak,
        cooling_rate: cooling.abs(),
        distance,
    }
}

/// Compute the thermal history at a fixed point `r = (0, d, 0)` as
/// the source moves past it.
pub fn thermal_history(
    src: &MovingPointSource,
    distance: f64,
    num_samples: usize,
) -> ThermalHistory {
    let v = src.scan_speed;
    let span = 10.0 * distance.max(1.0e-6) / v.max(1.0e-9);
    let dt = 2.0 * span / (num_samples as f64 - 1.0).max(1.0);
    let mut times = Vec::with_capacity(num_samples);
    let mut temps = Vec::with_capacity(num_samples);
    let mut peak_temp = f64::NEG_INFINITY;
    let mut peak_time = 0.0;
    for i in 0..num_samples {
        let t = -span + i as f64 * dt;
        let xi = -v * t;
        let temp = rosenthal_2d(src, xi, distance, 0.0);
        times.push(t);
        temps.push(temp);
        if temp > peak_temp {
            peak_temp = temp;
            peak_time = t;
        }
    }
    // Finite-difference cooling rate at the peak.
    let idx_peak = times
        .iter()
        .position(|&t| (t - peak_time).abs() < 1.0e-12)
        .unwrap_or(0);
    let cool = if idx_peak > 0 && idx_peak < temps.len() - 1 {
        let dt_local = times[idx_peak + 1] - times[idx_peak - 1];
        (temps[idx_peak + 1] - temps[idx_peak - 1]) / dt_local
    } else {
        0.0
    };
    ThermalHistory {
        times,
        temperatures: temps,
        peak_temperature: peak_temp,
        time_at_peak: peak_time,
        cooling_rate_at_peak: cool.abs(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn steel_point() -> MovingPointSource {
        MovingPointSource {
            linear_heat_input: 200.0 / 0.005, // P=200W, v=5 mm/s
            scan_speed: 0.005,
            thermal_conductivity: 30.0,
            thermal_diffusivity: 1.0e-5,
            preheat_temperature: 298.0,
        }
    }

    #[test]
    fn rosenthal_peak_above_preheat() {
        let s = steel_point();
        let t = rosenthal_2d(&s, 0.0, 0.005, 0.0);
        assert!(t > s.preheat_temperature);
    }

    #[test]
    fn rosenthal_drops_with_distance() {
        let s = steel_point();
        let t1 = rosenthal_2d(&s, 0.0, 0.001, 0.0);
        let t2 = rosenthal_2d(&s, 0.0, 0.01, 0.0);
        assert!(t1 > t2);
    }

    #[test]
    fn cooling_rate_positive() {
        let s = steel_point();
        let c = cooling_rate(&s, 0.005);
        assert!(c.cooling_rate > 0.0);
        assert!(c.peak_temperature > s.preheat_temperature);
    }

    #[test]
    fn thermal_history_peak_finite() {
        let s = steel_point();
        let h = thermal_history(&s, 0.005, 51);
        assert!(h.peak_temperature > s.preheat_temperature);
        assert!(h.cooling_rate_at_peak > 0.0);
        assert_eq!(h.times.len(), 51);
    }
}
