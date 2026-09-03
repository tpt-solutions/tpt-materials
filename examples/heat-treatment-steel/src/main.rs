//! Heat-treatment example.
//!
//! Simulates a 4-step heat-treatment schedule for an
//! AISI 4140-style steel (≈0.4 wt% C, M_s ≈ 600 K):
//!
//! 1. **Austenitise** + anneal (slow cool to room temp).
//! 2. **Austenitise** + oil quench to room temp.
//! 3. **Temper** the quenched structure for 1 hour at 600 °C.
//! 4. **Age** the tempered structure (peak-aged hardness
//!    increment).
//!
//! Each step is fed to [`simulate`] and prints the resulting
//! phase fractions and Vickers hardness.

use tpt_mat_heat::{
    hardness_from_fractions, hardness_martensite, simulate, AgingParams,
    AnnealingParams, HeatTreatmentProcess, PhaseHardness, QuenchingParams,
    TemperingParams,
};
use tpt_mat_heat_treatment as tpt_mat_heat;

fn report(label: &str, r: &tpt_mat_heat::HeatTreatmentResult) {
    println!(
        "{label:<9} | f_M={:>5.3} f_B={:>5.3} f_P={:>5.3} f_F={:>5.3} f_A={:>5.3} | HV = {:>6.1}",
        r.martensite_fraction,
        r.bainite_fraction,
        r.pearlite_fraction,
        r.ferrite_fraction,
        r.austenite_fraction,
        r.hardness_hv,
    );
}

fn main() {
    println!("AISI 4140-style heat-treatment schedule (0.4 wt% C, M_s = 600 K)");
    println!("----------+---------------------------------------------------+----------");
    println!("Step      | Phase fractions                                   | Hardness");
    println!("          | (M,B,P,F,A)                                       | (HV)");
    println!("----------+---------------------------------------------------+----------");

    let anneal = HeatTreatmentProcess::Annealing(AnnealingParams {
        austenitise_temperature: 1100.0,
        carbon_wt_pct: 0.4,
        equilibrium_ferrite_fraction: 0.5,
    });
    let r_anneal = simulate(&anneal).unwrap();
    report("Anneal", &r_anneal);

    let quench = HeatTreatmentProcess::Quenching(QuenchingParams {
        austenitise_temperature: 1100.0,
        quench_temperature: 300.0,
        ms_temperature: 600.0,
        km_alpha: 0.011,
        carbon_wt_pct: 0.4,
    });
    let r_quench = simulate(&quench).unwrap();
    report("Quench", &r_quench);

    let temper = HeatTreatmentProcess::Tempering(TemperingParams {
        initial_martensite: r_quench.martensite_fraction,
        temper_temperature: 873.0,
        hold_time: 3600.0,
        tempering_rate: 1.0e-3,
    });
    let r_temper = simulate(&temper).unwrap();
    report("Temper", &r_temper);

    let age = HeatTreatmentProcess::Aging(AgingParams {
        aging_temperature: 723.0,
        hold_time: 3600.0,
        peak_time: 3600.0,
        peak_age_hardening: 80.0,
    });
    let r_age = simulate(&age).unwrap();
    report("Age", &r_age);

    println!("----------+---------------------------------------------------+----------");

    // Cross-check the rule-of-mixtures and Krauss regression.
    let table = PhaseHardness::default();
    let hv_check = hardness_from_fractions(
        &[r_quench.martensite_fraction, 0.0, 0.0, 0.0, r_quench.austenite_fraction],
        &table,
        0.0,
    );
    println!(
        "Rule-of-mixtures HV (quench) = {:.1} (martensite-only Krauss = {:.1})",
        hv_check,
        hardness_martensite(0.4)
    );
}