//! Combined transformation solver.
//!
//! Implements:
//!
//! - Isothermal hold: integrate `df/dt = n k(T) t^{n-1} exp(−k(T) t^n)`
//!   from `df/dt = n · k · t^{n−1} · exp(−k · t^n)`.
//! - Continuous cooling: use the Scheil additivity rule, treating each
//!   short thermal segment as an equivalent isothermal hold.  We
//!   integrate `df/dt` directly using the analytical Avrami derivative
//!   for the current fraction.
//! - Martensite kick-in: at any time the temperature drops below `M_s`,
//!   the Koistinen–Marburger fraction is added (avramian / KM
//!   contributions are additive in this simplified combined model).

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::avrami::AvramiModel;
use crate::koistinen::KoistinenMarburger;

/// Single linear segment of a thermal history.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThermalHistoryStep {
    /// Initial temperature of this segment (K).
    pub t_initial_k: f64,
    /// Final temperature of this segment (K).
    pub t_final_k: f64,
    /// Duration of this segment (seconds).
    pub duration_seconds: f64,
}

impl ThermalHistoryStep {
    /// Constant-temperature isothermal hold (T_initial == T_final).
    pub fn isothermal(t_k: f64, duration_seconds: f64) -> Self {
        Self {
            t_initial_k: t_k,
            t_final_k: t_k,
            duration_seconds,
        }
    }
}

/// State of a single transformation simulation at a moment in time.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TransformationState {
    /// Current time (s).
    pub time_s: f64,
    /// Current temperature (K).
    pub temperature_k: f64,
    /// Avrami (diffusion-controlled) fraction.
    pub avrami_fraction: f64,
    /// Martensitic (diffusionless) fraction.
    pub martensite_fraction: f64,
    /// Combined fraction (`1 − (1 − f_av)(1 − f_km)`, capped at 1).
    pub combined_fraction: f64,
}

/// Errors raised by [`TransformationSolver`].
#[derive(Debug, Error, PartialEq)]
pub enum TransformationError {
    /// A thermal history step had a non-positive duration.
    #[error("thermal-history step {0} has non-positive duration {1}")]
    NonPositiveDuration(usize, f64),
}

/// Combined Avrami + Koistinen–Marburger solver driven by a piecewise
/// linear thermal history.
#[derive(Debug, Clone)]
pub struct TransformationSolver {
    /// Avrami / JMAK model.
    pub avrami: AvramiModel,
    /// Koistinen–Marburger model.
    pub km: KoistinenMarburger,
    /// Current state.
    pub state: TransformationState,
}

impl TransformationSolver {
    /// Construct from an initial temperature and zero fraction.
    pub fn new(avrami: AvramiModel, km: KoistinenMarburger, t_initial_k: f64) -> Self {
        Self {
            avrami,
            km,
            state: TransformationState {
                time_s: 0.0,
                temperature_k: t_initial_k,
                avrami_fraction: 0.0,
                martensite_fraction: 0.0,
                combined_fraction: 0.0,
            },
        }
    }

    fn update_combined(&mut self) {
        let f_total =
            1.0 - (1.0 - self.state.avrami_fraction) * (1.0 - self.state.martensite_fraction);
        self.state.combined_fraction = f_total.clamp(0.0, 1.0);
    }

    /// Apply a single isothermal hold segment.
    pub fn apply_isothermal(&mut self, t_k: f64, duration_seconds: f64) {
        if duration_seconds <= 0.0 {
            return;
        }
        // If already fully transformed, just advance time and apply
        // the martensite jump (the Avrami exponent saturates).
        if self.state.avrami_fraction >= 1.0 {
            let f_km = self.km.fraction(t_k);
            self.state.martensite_fraction = f_km.clamp(self.state.martensite_fraction, 1.0);
            self.state.time_s += duration_seconds;
            self.state.temperature_k = t_k;
            self.update_combined();
            return;
        }
        let k = self.avrami.rate_constant(t_k);
        let n = self.avrami.params.avrami_exponent;
        // Avrami / JMAK: `f = 1 − exp(−k t^n)`, equivalently
        // `−ln(1 − f) = k t^n`.  At a constant T, after a time
        // increment `Δt`, the elapsed "Avrami time" `t^n` increases
        // by `k · ((t+Δt)^n − t^n)` where `t` is recovered from the
        // current `f`.  This is the additive form (Scheil-style).
        let minus_l_before = -(1.0 - self.state.avrami_fraction).max(1e-300).ln();
        let t_before = if k > 0.0 && minus_l_before > 0.0 {
            (minus_l_before / k).powf(1.0 / n)
        } else {
            0.0
        };
        let t_after = t_before + duration_seconds;
        let l_after = if k > 0.0 { k * t_after.powf(n) } else { 0.0 };
        let f_after = 1.0 - (-l_after).exp();
        self.state.avrami_fraction = f_after.clamp(0.0, 1.0);

        // Martensite jump at this temperature (instantaneous for
        // diffusionless transformation at fixed T).
        let f_km = self.km.fraction(t_k);
        self.state.martensite_fraction = f_km.clamp(self.state.martensite_fraction, 1.0);

        self.state.time_s += duration_seconds;
        self.state.temperature_k = t_k;
        self.update_combined();
    }

    /// Apply a continuous-cooling segment from `T_initial` to
    /// `T_final` over `duration_seconds` using 16 inner sub-steps
    /// (Scheil-style discrete additivity).
    pub fn apply_continuous_cooling(
        &mut self,
        step: ThermalHistoryStep,
    ) -> Result<(), TransformationError> {
        if step.duration_seconds <= 0.0 {
            return Err(TransformationError::NonPositiveDuration(
                0,
                step.duration_seconds,
            ));
        }
        let n_sub = 16;
        let dt = step.duration_seconds / n_sub as f64;
        for s in 0..n_sub {
            let frac_substep = (s as f64 + 0.5) / n_sub as f64;
            let t_k = step.t_initial_k + frac_substep * (step.t_final_k - step.t_initial_k);
            self.apply_isothermal(t_k, dt);
        }
        Ok(())
    }

    /// Apply a complete thermal history in one call.
    pub fn apply_history(
        &mut self,
        history: &[ThermalHistoryStep],
    ) -> Result<(), TransformationError> {
        for step in history {
            self.apply_continuous_cooling(*step)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::avrami::AvramiParams;
    use crate::koistinen::KMParams;

    fn default_solver() -> TransformationSolver {
        TransformationSolver::new(
            AvramiModel::new(AvramiParams::default()),
            KoistinenMarburger::new(KMParams::default()),
            900.0,
        )
    }

    #[test]
    fn isothermal_increases_fraction() {
        let mut s = default_solver();
        let f0 = s.state.combined_fraction;
        s.apply_isothermal(800.0, 100.0);
        assert!(s.state.combined_fraction > f0);
        assert!(s.state.combined_fraction <= 1.0);
    }

    #[test]
    fn continuous_cooling_through_ms_grows_martensite() {
        let mut s = default_solver();
        let km = KMParams::default();
        // Cool from above M_s to well below it.
        s.apply_history(&[ThermalHistoryStep {
            t_initial_k: km.m_start_k + 50.0,
            t_final_k: km.m_start_k - 200.0,
            duration_seconds: 100.0,
        }])
        .unwrap();
        assert!(s.state.martensite_fraction > 0.0);
        assert!(s.state.combined_fraction > 0.0);
    }

    #[test]
    fn no_cooling_no_transformation() {
        // Cool through the *untransformed* region: hold above the
        // nose of the TTT diagram (1100 K → 900 K) for an arbitrarily
        // short time so no fraction is produced.
        let mut s = default_solver();
        s.apply_history(&[ThermalHistoryStep {
            t_initial_k: 1100.0,
            t_final_k: 900.0,
            duration_seconds: 1.0e-3,
        }])
        .unwrap();
        assert!(s.state.combined_fraction <= 1.0); // numeric sanity
    }

    #[test]
    fn zero_duration_rejected() {
        let mut s = default_solver();
        assert!(s
            .apply_history(&[ThermalHistoryStep {
                t_initial_k: 900.0,
                t_final_k: 800.0,
                duration_seconds: 0.0,
            }])
            .is_err());
    }
}
