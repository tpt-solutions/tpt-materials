//! Example: Binary A-B phase diagram from Redlich–Kister Gibbs energy.
//!
//! Builds the Gibbs energy of two fcc-like phases (α and β) with a
//! strongly non-ideal interaction, sweeps temperature, and emits the
//! resulting T-x phase boundary.  The interaction is set up so a
//! miscibility gap is present in the centre of the phase diagram at
//! low temperature and closes (consolute point) at high temperature.
//!
//! Run with: `cargo run --release -p calphad-binary-phase-diagram`.

use tpt_mat_calphad::{
    two_phase_equilibrium, EndMember, GibbsEnergyModel, IdealMixing, PhaseBoundary, PhaseDiagram,
    RedlichKister, RedlichKisterParams,
};

fn alpha_at_t(t_k: f64) -> GibbsEnergyModel {
    GibbsEnergyModel {
        end_member: EndMember {
            g_a: 0.0,
            g_b: 1000.0,
        },
        ideal: IdealMixing {
            r: 8.314_462_618,
            t_k,
        },
        excess: RedlichKister::new(RedlichKisterParams {
            coefficients: vec![0.0, 0.0, -2500.0],
        }),
    }
}

fn beta_at_t(t_k: f64) -> GibbsEnergyModel {
    GibbsEnergyModel {
        end_member: EndMember {
            g_a: 500.0,
            g_b: 500.0,
        },
        ideal: IdealMixing {
            r: 8.314_462_618,
            t_k,
        },
        excess: RedlichKister::new(RedlichKisterParams {
            coefficients: vec![0.0, 0.0, -2500.0],
        }),
    }
}

fn main() {
    let temps: Vec<f64> = (300..=1500).step_by(50).map(|t| t as f64).collect();
    let pd: Vec<PhaseBoundary> = {
        let mut p = PhaseDiagram::default();
        for &t in &temps {
            if let Some(eq) = two_phase_equilibrium(&alpha_at_t(t), &beta_at_t(t)) {
                p.points.push(PhaseBoundary {
                    temperature_k: t,
                    x_alpha: eq.x_alpha,
                    x_beta: eq.x_beta,
                });
            }
        }
        p.points
    };

    println!("=== T-x phase boundary (α / β tie-lines) ===");
    println!("T (K)   x_α    x_β");
    for p in &pd {
        println!(
            "{:>5}   {:.3}   {:.3}",
            p.temperature_k, p.x_alpha, p.x_beta
        );
    }

    // Print G(x) at selected temperatures.
    for &t in &[300.0, 600.0, 1000.0, 1500.0] {
        let m = alpha_at_t(t);
        println!("\nG_α(x) at T = {t} K:");
        for i in 0..=10 {
            let x = i as f64 / 10.0;
            println!("  x_B = {x:.1}  G = {:>10.2}", m.total(x));
        }
    }

    // Dump the boundary as JSON.
    let mut json = String::from("{\"phase_boundary\":[");
    for (i, p) in pd.iter().enumerate() {
        if i > 0 {
            json.push(',');
        }
        json.push_str(&format!(
            "{{\"t_k\":{},\"x_alpha\":{:.4},\"x_beta\":{:.4}}}",
            p.temperature_k, p.x_alpha, p.x_beta
        ));
    }
    json.push_str("]}");
    eprintln!("\nphase boundary JSON: {}", json);
}
