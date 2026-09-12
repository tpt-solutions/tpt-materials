use tpt_mat_crystal_plasticity::SymmetricFourthOrder;
use tpt_mat_rve::{HomogenizationScheme, Rve, RveGrain};
use tpt_math_linalg_fixed::Vec6;

fn isotropic_oriented(label: &str, volume_fraction: f64, e: f64, nu: f64) -> RveGrain {
    let c = SymmetricFourthOrder::isotropic(e, nu);
    RveGrain::new(
        label,
        volume_fraction,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        c,
    )
}

fn anisotropic_fcc(label: &str, volume_fraction: f64) -> RveGrain {
    let c = SymmetricFourthOrder::cubic(200_000.0, 120_000.0, 80_000.0);
    RveGrain::new(
        label,
        volume_fraction,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        c,
    )
}

#[test]
fn hill_mandel_voigt_end_to_end() {
    let rve = Rve::new(vec![
        anisotropic_fcc("a", 0.3),
        isotropic_oriented("b", 0.7, 180_000.0, 0.28),
    ]);
    let c_avg = rve.homogenize(HomogenizationScheme::Voigt);
    let e = Vec6::new(1.0e-3, -0.5e-3, 0.0, 0.0, 0.0, 0.5e-3);
    let sigma_avg = c_avg.contract(e);
    let sigma_a = rve.grains[0].stiffness.contract(e);
    let sigma_b = rve.grains[1].stiffness.contract(e);
    let macro_energy = sigma_avg.double_dot(e);
    let micro_energy = 0.3 * sigma_a.double_dot(e) + 0.7 * sigma_b.double_dot(e);
    let rel = (macro_energy - micro_energy).abs() / macro_energy.abs().max(1e-30);
    assert!(rel < 1.0e-12, "Voigt Hill–Mandel mismatch: {rel:.3e}");
}

#[test]
fn hill_mandel_reuss_end_to_end() {
    let c_a = SymmetricFourthOrder::cubic(200_000.0, 120_000.0, 80_000.0);
    let c_b = SymmetricFourthOrder::isotropic(180_000.0, 0.28);
    let phases = vec![(c_a.clone(), 0.3), (c_b.clone(), 0.7)];
    let s_avg = tpt_mat_rve::homogenization::reuss(&phases);
    let s_a = c_a.compliance();
    let s_b = c_b.compliance();
    let sigma = Vec6::new(100.0, -50.0, 0.0, 0.0, 0.0, 25.0);
    let apply = |s: [[f64; 6]; 6]| {
        let mut out = [0.0_f64; 6];
        for i in 0..6 {
            let mut acc = 0.0;
            for j in 0..6 {
                acc += s[i][j] * sigma.data[j];
            }
            out[i] = acc;
        }
        Vec6::new(out[0], out[1], out[2], out[3], out[4], out[5])
    };
    let e_avg = apply(s_avg);
    let e_a = apply(s_a).scale(0.3);
    let e_b = apply(s_b).scale(0.7);
    let macro_energy = sigma.double_dot(e_avg);
    let micro_energy = sigma.double_dot(e_a + e_b);
    let rel = (macro_energy - micro_energy).abs() / macro_energy.abs().max(1e-30);
    assert!(rel < 1.0e-9, "Reuss Hill–Mandel mismatch: {rel:.3e}");
}

#[test]
fn hill_mandel_self_consistent_end_to_end() {
    let c_a = SymmetricFourthOrder::cubic(200_000.0, 120_000.0, 80_000.0);
    let c_b = SymmetricFourthOrder::isotropic(180_000.0, 0.28);
    let rve = Rve::new(vec![
        RveGrain::new(
            "a",
            0.3,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            c_a.clone(),
        ),
        RveGrain::new(
            "b",
            0.7,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            c_b.clone(),
        ),
    ]);
    let c_sc = rve.homogenize(HomogenizationScheme::SelfConsistent);
    let e = Vec6::new(1.0e-3, -0.5e-3, 0.0, 0.0, 0.0, 0.5e-3);
    let sigma_avg = c_sc.contract(e);
    let sigma_a = c_a.contract(e);
    let sigma_b = c_b.contract(e);
    let macro_energy = sigma_avg.double_dot(e);
    let micro_energy = 0.3 * sigma_a.double_dot(e) + 0.7 * sigma_b.double_dot(e);
    let rel = (macro_energy - micro_energy).abs() / macro_energy.abs().max(1e-30);
    assert!(rel < 1.0, "Self-consistent Hill–Mandel mismatch: {rel:.3e}");
}
