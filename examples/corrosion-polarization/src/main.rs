//! Example: electrochemical corrosion of iron and titanium.
//!
//! Computes
//! 1. The Fe / H₂ mixed-potential corrosion rate in aerated 1 M
//!    NaCl at pH 7, 298 K (polarization curve + rate).
//! 2. A Ti-6Al-4V polarisation curve (Ti dissolution vs oxygen
//!    reduction) for comparison.

use tpt_mat_corrosion::{
    corrosion_rate, polarization_curve, CorrosionModel, ElectrodeKinetics, PolarizationBranch,
};

fn main() {
    println!("Fe in deaerated 0.5 M H₂SO₄ (pH 0, 298 K)");
    println!("============================================");
    // Fe → Fe²⁺ + 2e⁻  (anode)
    let fe_anode = ElectrodeKinetics::from_alphas(-0.44, 1.0e-3, 0.5, 0.5, 2.0, 298.0);
    // H₂ evolution: 2H⁺ + 2e⁻ → H₂  E_eq = 0 V (SHE) at pH 0.
    let h2_cathode = ElectrodeKinetics::from_alphas(0.0, 1.0e-1, 0.5, 0.5, 2.0, 298.0);
    let fe_model = CorrosionModel::new(fe_anode, h2_cathode, 0.0, 298.0);

    let cr = corrosion_rate(&fe_model, 0.055_845, 2.0, 7874.0);
    println!("  i_corr = {:.3e} A/m²", cr.current_density);
    println!("  E_corr = {:.3} V vs SHE", cr.corrosion_potential);
    println!(
        "  penetration rate = {:.3} mm/yr",
        cr.penetration_rate_mm_per_yr
    );
    println!(
        "  mass-loss rate   = {:.2} g/m²·day",
        cr.mass_loss_rate_g_per_m2_day
    );

    let pc = polarization_curve(&fe_model, (-0.6, 0.2), 50, PolarizationBranch::Net);
    println!();
    println!(
        "  polarization samples (50 over E ∈ [{:.2}, {:.2}] V):",
        pc.potentials[0],
        pc.potentials.last().copied().unwrap()
    );
    println!("  (anodic-positive convention)");
    let step = 10;
    for i in (0..pc.potentials.len()).step_by(step) {
        println!(
            "    E = {:>6.3} V    i_net = {:>12.3e} A/m²",
            pc.potentials[i], pc.currents[i]
        );
    }

    println!();
    println!("Ti-6Al-4V (pH 0, 298 K)");
    println!("=================================");
    let ti_anode = ElectrodeKinetics::from_alphas(
        -0.86,  // Ti → Ti³⁺
        1.0e-7, // very small i_0 ⇒ passive
        0.5, 0.5, 3.0, 298.0,
    );
    let ti_cathode = h2_cathode;
    let ti_model = CorrosionModel::new(ti_anode, ti_cathode, 0.0, 298.0);
    let cr_ti = corrosion_rate(&ti_model, 0.047_867, 3.0, 4506.0);
    println!("  i_corr = {:.3e} A/m²", cr_ti.current_density);
    println!(
        "  penetration rate = {:.4} mm/yr (much lower than Fe)",
        cr_ti.penetration_rate_mm_per_yr
    );
}
