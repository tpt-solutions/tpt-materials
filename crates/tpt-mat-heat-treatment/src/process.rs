//! Heat-treatment process definitions and the [`simulate`] driver.

use serde::{Deserialize, Serialize};
use thiserror::Error;
use tpt_mat_phase_transform::{AvramiModel, AvramiParams, KMParams, KoistinenMarburger};

use crate::hardness::{hardness_from_fractions, PhaseHardness};

/// Heat-treatment process variants.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum HeatTreatmentProcess {
    /// Annealing — slow cool from austenitising temperature.
    Annealing(AnnealingParams),
    /// Quenching — rapid cool to a quench medium temperature.
    Quenching(QuenchingParams),
    /// Tempering — sub-critical reheat after a quench.
    Tempering(TemperingParams),
    /// Aging — low-temperature precipitation hardening.
    Aging(AgingParams),
}

/// Parameters for an annealing treatment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AnnealingParams {
    /// Austenitising temperature (K).
    pub austenitise_temperature: f64,
    /// Carbon content (wt%) of the alloy.
    pub carbon_wt_pct: f64,
    /// Equilibrium ferrite fraction at room temperature
    /// (`0 ≤ f ≤ 1`).
    pub equilibrium_ferrite_fraction: f64,
}

/// Parameters for a quenching treatment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct QuenchingParams {
    /// Austenitising temperature (K).
    pub austenitise_temperature: f64,
    /// Quench-medium temperature (K).
    pub quench_temperature: f64,
    /// Martensite-start temperature (K).
    pub ms_temperature: f64,
    /// Koistinen–Marburger `α` parameter (typical `0.011` for steels).
    pub km_alpha: f64,
    /// Carbon content (wt%) — drives martensite hardness.
    pub carbon_wt_pct: f64,
}

/// Parameters for a tempering treatment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TemperingParams {
    /// Initial martensite fraction after the prior quench (`0..=1`).
    pub initial_martensite: f64,
    /// Tempering temperature (K).
    pub temper_temperature: f64,
    /// Hold time at temperature (s).
    pub hold_time: f64,
    /// Tempering rate constant `k_t` (1/s).
    pub tempering_rate: f64,
}

/// Parameters for an aging treatment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AgingParams {
    /// Aging temperature (K).
    pub aging_temperature: f64,
    /// Hold time (s).
    pub hold_time: f64,
    /// Time to peak hardness `t_peak` (s).
    pub peak_time: f64,
    /// Peak age-hardening contribution `ΔHV_peak` (HV).
    pub peak_age_hardening: f64,
}

/// Result of a heat-treatment simulation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HeatTreatmentResult {
    /// Process variant that produced this result.
    pub process: HeatTreatmentProcess,
    /// Final martensite fraction.
    pub martensite_fraction: f64,
    /// Final bainite fraction.
    pub bainite_fraction: f64,
    /// Final pearlite fraction.
    pub pearlite_fraction: f64,
    /// Final ferrite fraction.
    pub ferrite_fraction: f64,
    /// Final retained-austenite fraction.
    pub austenite_fraction: f64,
    /// Predicted Vickers hardness (HV).
    pub hardness_hv: f64,
    /// Martensite hardness contribution (HV) used in the mix,
    /// if applicable.
    pub martensite_hardness_hv: f64,
}

/// Errors produced by the heat-treatment simulator.
#[derive(Debug, Error, PartialEq)]
pub enum HeatTreatmentError {
    /// Inputs violated a `0 ≤ x ≤ 1` constraint.
    #[error("fraction out of [0,1]: {0}")]
    InvalidFraction(String),
    /// Temperature is non-positive.
    #[error("non-physical temperature: {0}")]
    InvalidTemperature(String),
}

/// Simulate a heat-treatment process.
pub fn simulate(process: &HeatTreatmentProcess) -> Result<HeatTreatmentResult, HeatTreatmentError> {
    match process {
        HeatTreatmentProcess::Annealing(p) => simulate_annealing(p),
        HeatTreatmentProcess::Quenching(p) => simulate_quench(p),
        HeatTreatmentProcess::Tempering(p) => simulate_tempering(p),
        HeatTreatmentProcess::Aging(p) => simulate_aging(p),
    }
}

fn simulate_annealing(p: &AnnealingParams) -> Result<HeatTreatmentResult, HeatTreatmentError> {
    if p.austenitise_temperature <= 0.0 {
        return Err(HeatTreatmentError::InvalidTemperature(format!(
            "austenitise_temperature = {}",
            p.austenitise_temperature
        )));
    }
    if !(0.0..=1.0).contains(&p.equilibrium_ferrite_fraction) {
        return Err(HeatTreatmentError::InvalidFraction(format!(
            "equilibrium_ferrite_fraction = {}",
            p.equilibrium_ferrite_fraction
        )));
    }
    // Slow cool → equilibrium phase fractions:
    //   ferrite = equilibrium_ferrite_fraction
    //   pearlite = 1 - ferrite (residual carbon-bearing phases)
    //   martensite = 0
    let ferrite = p.equilibrium_ferrite_fraction;
    let pearlite = 1.0 - ferrite;
    let fractions = [0.0, 0.0, pearlite, ferrite, 0.0];
    let table = PhaseHardness::default();
    let hv = hardness_from_fractions(&fractions, &table, 0.0);
    Ok(HeatTreatmentResult {
        process: HeatTreatmentProcess::Annealing(*p),
        martensite_fraction: 0.0,
        bainite_fraction: 0.0,
        pearlite_fraction: pearlite,
        ferrite_fraction: ferrite,
        austenite_fraction: 0.0,
        hardness_hv: hv,
        martensite_hardness_hv: 0.0,
    })
}

fn simulate_quench(p: &QuenchingParams) -> Result<HeatTreatmentResult, HeatTreatmentError> {
    if p.ms_temperature <= 0.0 {
        return Err(HeatTreatmentError::InvalidTemperature(format!(
            "ms_temperature = {}",
            p.ms_temperature
        )));
    }
    // Koistinen–Marburger for the athermal martensite formation.
    let km = KoistinenMarburger::new(KMParams {
        m_start_k: p.ms_temperature,
        alpha_per_k: p.km_alpha,
        retained_at_m_start: 0.01,
    });
    let martensite = km.fraction(p.quench_temperature);
    let retained_austenite = (1.0_f64 - martensite).max(0.0);
    let fractions = [martensite, 0.0, 0.0, 0.0, retained_austenite];
    let martensite_hv = crate::hardness::hardness_martensite(p.carbon_wt_pct);
    // Custom phase hardness table that scales martensite hardness
    // with carbon content.
    let table = PhaseHardness {
        martensite: martensite_hv,
        ..PhaseHardness::default()
    };
    let hv = hardness_from_fractions(&fractions, &table, 0.0);
    Ok(HeatTreatmentResult {
        process: HeatTreatmentProcess::Quenching(*p),
        martensite_fraction: martensite,
        bainite_fraction: 0.0,
        pearlite_fraction: 0.0,
        ferrite_fraction: 0.0,
        austenite_fraction: retained_austenite,
        hardness_hv: hv,
        martensite_hardness_hv: martensite_hv,
    })
}

fn simulate_tempering(p: &TemperingParams) -> Result<HeatTreatmentResult, HeatTreatmentError> {
    if !(0.0..=1.0).contains(&p.initial_martensite) {
        return Err(HeatTreatmentError::InvalidFraction(format!(
            "initial_martensite = {}",
            p.initial_martensite
        )));
    }
    if p.temper_temperature <= 0.0 {
        return Err(HeatTreatmentError::InvalidTemperature(format!(
            "temper_temperature = {}",
            p.temper_temperature
        )));
    }
    // First-order decay: martensite -> tempered martensite.
    // We model the tempered fraction as carbide + ferrite, treated
    // as ferrite-equivalent for hardness purposes.
    let tempered_fraction = p.initial_martensite * (1.0 - (-p.tempering_rate * p.hold_time).exp());
    let martensite = p.initial_martensite - tempered_fraction;
    let ferrite = tempered_fraction; // tempered martensite ≈ ferrite + carbides
    let fractions = [martensite, 0.0, 0.0, ferrite, 0.0];
    let table = PhaseHardness::default();
    let hv = hardness_from_fractions(&fractions, &table, 0.0);
    Ok(HeatTreatmentResult {
        process: HeatTreatmentProcess::Tempering(*p),
        martensite_fraction: martensite,
        bainite_fraction: 0.0,
        pearlite_fraction: 0.0,
        ferrite_fraction: ferrite,
        austenite_fraction: 0.0,
        hardness_hv: hv,
        martensite_hardness_hv: 0.0,
    })
}

fn simulate_aging(p: &AgingParams) -> Result<HeatTreatmentResult, HeatTreatmentError> {
    if p.aging_temperature <= 0.0 {
        return Err(HeatTreatmentError::InvalidTemperature(format!(
            "aging_temperature = {}",
            p.aging_temperature
        )));
    }
    // Gaussian age-hardening curve centred at peak_time.
    // ΔHV(t) = ΔHV_peak exp(−½ ((t − t_peak)/τ)²)
    // with τ = peak_time / 2 (so the curve decays by peak/2).
    let tau = p.peak_time.max(1.0) / 2.0;
    let z = (p.hold_time - p.peak_time) / tau;
    let extra_hv = p.peak_age_hardening * (-0.5 * z * z).exp();

    // The matrix is a fixed starting state (martensite-as-quenched,
    // then tempered, then aged): we use a 100% martensite proxy to
    // demonstrate the age-hardening bump on a high-strength substrate.
    let fractions = [1.0, 0.0, 0.0, 0.0, 0.0];
    let table = PhaseHardness::default();
    let hv = hardness_from_fractions(&fractions, &table, extra_hv);
    Ok(HeatTreatmentResult {
        process: HeatTreatmentProcess::Aging(*p),
        martensite_fraction: 1.0,
        bainite_fraction: 0.0,
        pearlite_fraction: 0.0,
        ferrite_fraction: 0.0,
        austenite_fraction: 0.0,
        hardness_hv: hv,
        martensite_hardness_hv: table.martensite,
    })
}

// `AvramiModel` / `AvramiParams` are not directly used here —
// re-export to suppress dead-code warnings if the crate later
// builds on these symbols.
#[allow(dead_code)]
fn _avrami_link(_a: &AvramiModel, _p: &AvramiParams) {}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn annealing_yields_equilibrium_fractions() {
        let p = AnnealingParams {
            austenitise_temperature: 1100.0,
            carbon_wt_pct: 0.4,
            equilibrium_ferrite_fraction: 0.5,
        };
        let r = simulate(&HeatTreatmentProcess::Annealing(p)).unwrap();
        assert!(approx(r.ferrite_fraction, 0.5, 1.0e-12));
        assert!(approx(r.pearlite_fraction, 0.5, 1.0e-12));
        assert!(approx(r.martensite_fraction, 0.0, 1.0e-12));
    }

    #[test]
    fn quenching_below_ms_yields_martensite() {
        // Quench to 200 K, M_s = 600 K, α = 0.011:
        // f = 1 − exp(−α (M_s − T)) = 1 − exp(−0.011 * 400) = 1 − exp(−4.4) ≈ 0.9877.
        let p = QuenchingParams {
            austenitise_temperature: 1100.0,
            quench_temperature: 200.0,
            ms_temperature: 600.0,
            km_alpha: 0.011,
            carbon_wt_pct: 0.4,
        };
        let r = simulate(&HeatTreatmentProcess::Quenching(p)).unwrap();
        let expected = 1.0 - (-0.011_f64 * 400.0).exp();
        assert!(
            approx(r.martensite_fraction, expected, 1.0e-9),
            "f_M = {} (expected {})",
            r.martensite_fraction,
            expected
        );
        assert!(r.hardness_hv > 400.0);
    }

    #[test]
    fn quenching_above_ms_retains_austenite() {
        // Quench to 700 K, M_s = 600 K — no martensite forms.
        let p = QuenchingParams {
            austenitise_temperature: 1100.0,
            quench_temperature: 700.0,
            ms_temperature: 600.0,
            km_alpha: 0.011,
            carbon_wt_pct: 0.4,
        };
        let r = simulate(&HeatTreatmentProcess::Quenching(p)).unwrap();
        assert!(approx(r.martensite_fraction, 0.0, 1.0e-12));
        assert!(approx(r.austenite_fraction, 1.0, 1.0e-12));
    }

    #[test]
    fn tempering_reduces_martensite_exponentially() {
        let p = TemperingParams {
            initial_martensite: 1.0,
            temper_temperature: 873.0,
            hold_time: 3600.0,
            tempering_rate: 1.0e-3,
        };
        let r = simulate(&HeatTreatmentProcess::Tempering(p)).unwrap();
        let tempered = 1.0 - (-1.0e-3_f64 * 3600.0).exp();
        assert!(
            approx(r.ferrite_fraction, tempered, 1.0e-9),
            "ferrite = {} (expected {})",
            r.ferrite_fraction,
            tempered
        );
        assert!(r.martensite_fraction < r.ferrite_fraction);
    }

    #[test]
    fn aging_produces_peak_hardness_at_peak_time() {
        let p_at_peak = AgingParams {
            aging_temperature: 450.0,
            hold_time: 3600.0,
            peak_time: 3600.0,
            peak_age_hardening: 100.0,
        };
        let p_off = AgingParams {
            hold_time: 0.0,
            ..p_at_peak
        };
        let r_peak = simulate(&HeatTreatmentProcess::Aging(p_at_peak)).unwrap();
        let r_off = simulate(&HeatTreatmentProcess::Aging(p_off)).unwrap();
        assert!(r_peak.hardness_hv > r_off.hardness_hv);
    }

    #[test]
    fn invalid_fraction_returns_error() {
        let p = AnnealingParams {
            austenitise_temperature: 1100.0,
            carbon_wt_pct: 0.4,
            equilibrium_ferrite_fraction: 1.5,
        };
        let r = simulate(&HeatTreatmentProcess::Annealing(p));
        assert!(matches!(r, Err(HeatTreatmentError::InvalidFraction(_))));
    }

    #[test]
    fn invalid_temperature_returns_error() {
        let p = QuenchingParams {
            austenitise_temperature: 1100.0,
            quench_temperature: 200.0,
            ms_temperature: -1.0,
            km_alpha: 0.011,
            carbon_wt_pct: 0.4,
        };
        let r = simulate(&HeatTreatmentProcess::Quenching(p));
        assert!(matches!(r, Err(HeatTreatmentError::InvalidTemperature(_))));
    }
}
