//! Example: Voigt / Reuss / Hashin–Shtrikman bounds for a
//! two-phase Al + steel composite.
//!
//! Demonstrates [`tpt_mat_homogenization::voigt_reuss_bounds`] and
//! [`tpt_mat_homogenization::hashin_shtrikman_two_phase`].
//!
//! Run with: `cargo run --release -p homogenization-voigt-reuss`.

use tpt_mat_crystal_plasticity::SymmetricFourthOrder;
use tpt_mat_homogenization::{
    g_from_e_nu, hashin_shtrikman_two_phase, k_from_e_nu, voigt_reuss_bounds,
};

fn main() {
    // Phase A: aluminum-like  (E = 70 GPa, ν = 0.33)
    // Phase B: steel-like     (E = 200 GPa, ν = 0.30)
    let e_a = 70_000.0_f64;
    let nu_a = 0.33_f64;
    let e_b = 200_000.0_f64;
    let nu_b = 0.30_f64;
    let k_a = k_from_e_nu(e_a, nu_a);
    let g_a = g_from_e_nu(e_a, nu_a);
    let k_b = k_from_e_nu(e_b, nu_b);
    let g_b = g_from_e_nu(e_b, nu_b);

    println!(
        "Phase A (Al-like):   E = {e_a:.0} MPa, ν = {nu_a}, K = {k_a:.1} MPa, G = {g_a:.1} MPa"
    );
    println!(
        "Phase B (steel-like): E = {e_b:.0} MPa, ν = {nu_b}, K = {k_b:.1} MPa, G = {g_b:.1} MPa"
    );

    // Sweep volume fraction of phase B.
    println!();
    println!(
        " f_B    | K_Voigt  K_Reuss  K_VRH   | G_Voigt  G_Reuss  G_VRH   | E_VRH (MPa)  ν_VRH"
    );
    println!(
        "--------+--------------------------+--------------------------+----------------+------"
    );
    for i in 0..=10 {
        let f_b = i as f64 / 10.0;
        let phases = [(k_a, g_a, 1.0 - f_b), (k_b, g_b, f_b)];
        let bounds = voigt_reuss_bounds(&phases);
        let hs = hashin_shtrikman_two_phase((e_a, nu_a, 1.0 - f_b), (e_b, nu_b, f_b));
        println!(
            " {:>5.2}  | {:>7.1}  {:>7.1}  {:>7.1} | {:>7.1}  {:>7.1}  {:>7.1} | {:>13.0}  {:>5.3}  [HS: K ∈ [{:.0}, {:.0}], G ∈ [{:.0}, {:.0}]]",
            f_b,
            bounds.k_voigt,
            bounds.k_reuss,
            bounds.k_vrh(),
            bounds.g_voigt,
            bounds.g_reuss,
            bounds.g_vrh(),
            bounds.e_vrh(),
            bounds.e_vrh() / (2.0 * bounds.g_vrh()) - 1.0,
            hs.k_lower, hs.k_upper, hs.g_lower, hs.g_upper,
        );
    }

    // Also compute the 6x6 Voigt average stiffness at f_B = 0.3.
    let c_a = SymmetricFourthOrder::isotropic(e_a, nu_a);
    let c_b = SymmetricFourthOrder::isotropic(e_b, nu_b);
    let f_b = 0.3_f64;
    let c_avg = tpt_mat_homogenization::voigt(&[(c_a.clone(), 1.0 - f_b), (c_b.clone(), f_b)]);
    println!();
    println!("Voigt-averaged 6x6 stiffness at f_B = {f_b}:");
    for i in 0..6 {
        let row: Vec<String> = (0..6)
            .map(|j| format!("{:>10.2}", c_avg.data[i][j]))
            .collect();
        println!("  [{}]", row.join(", "));
    }
}
