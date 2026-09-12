//! Regenerate the `test-data/golden/*` reference datasets.
//!
//! Every golden test in the workspace recomputes exactly what is written
//! here with a deterministic setup (fixed seeds/parameters), so running
//! `cargo run -p golden-generate` after a deliberate behavioural change
//! updates the committed baselines.  The `.gitignore` keeps any
//! `*.actual.json` comparison output out of the repository.

use std::fs;
use std::path::{Path, PathBuf};

use serde::Serialize;

use tpt_mat_battery::{capacity_fade_curve, default_active_material, BatteryChemistry};
use tpt_mat_corrosion::{
    corrosion_rate, polarization_curve, CorrosionModel, ElectrodeKinetics, PolarizationBranch,
};
use tpt_mat_crystal_plasticity::{
    solve_increment_single_point, CrystalPlasticityModel, SymmetricFourthOrder,
};
use tpt_mat_crystallography::{CrystalStructure, SlipSystem};
use tpt_mat_damage::{gtn_yield_function, update_porosity, GtnParams};
use tpt_mat_fatigue_micro::{predict_crack_initiation, FatigueCriterion};
use tpt_mat_hardening::{Hardening, HardeningState, VoceHardening, VoceParams};
use tpt_mat_phase_field::{BulkEnergy, PhaseFieldSolver, RegularSolutionParams};
use tpt_mat_rve::bishop_hill_taylor_factor_axis;
use tpt_mat_solidification::{AnisotropyMode, AnisotropyModel, SolidificationSolver};
use tpt_science::Grid2D;

/// Workspace `test-data/golden` directory (this example lives in
/// `examples/golden-generate`, two levels below the workspace root).
fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-data/golden")
}

fn write_json(name: &str, value: &impl Serialize) {
    let path = golden_dir().join(name);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).unwrap();
    }
    let text = serde_json::to_string_pretty(value).unwrap();
    fs::write(&path, format!("{text}\n")).unwrap();
    println!("wrote {}", path.display());
}

// ---------------------------------------------------------------------
// Deterministic pseudo-random helpers (shared with the golden tests).
// ---------------------------------------------------------------------

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    *seed as f64 / u64::MAX as f64
}

/// Uniform random unit vector (Marsaglia sphere sampling).
fn random_unit_axis(seed: &mut u64) -> [f64; 3] {
    loop {
        let u = 2.0 * lcg(seed) - 1.0;
        let v = 2.0 * lcg(seed) - 1.0;
        let s = u * u + v * v;
        if s < 1.0 && s > 0.0 {
            let f = 2.0 * (1.0 - s).sqrt();
            return [u * f, v * f, 1.0 - 2.0 * s];
        }
    }
}

// ---------------------------------------------------------------------
// 1. Taylor factors (Phase 2 golden data).
// ---------------------------------------------------------------------

#[derive(Serialize)]
struct AxisValue {
    name: String,
    axis: [f64; 3],
    m: f64,
}

#[derive(Serialize)]
struct RandomAxisValue {
    axis: [f64; 3],
    m: f64,
}

#[derive(Serialize)]
struct TaylorFile {
    label: String,
    slip_family: &'static str,
    n_random_total: usize,
    random_average: f64,
    axes: Vec<AxisValue>,
    random_axes: Vec<RandomAxisValue>,
}

fn taylor_dataset(
    label: &str,
    slip_family: &'static str,
    slips: Option<&[SlipSystem]>,
) -> TaylorFile {
    let named: [(&str, [f64; 3]); 6] = [
        ("001", [0.0, 0.0, 1.0]),
        ("110", [1.0, 1.0, 0.0]),
        ("111", [1.0, 1.0, 1.0]),
        ("012", [0.0, 1.0, 2.0]),
        ("112", [1.0, 1.0, 2.0]),
        ("123", [1.0, 2.0, 3.0]),
    ];
    let axes = named
        .iter()
        .map(|(name, axis)| AxisValue {
            name: name.to_string(),
            axis: *axis,
            m: bishop_hill_taylor_factor_axis(*axis, 1.0, slips).taylor_factor,
        })
        .collect();

    let n_total = 256usize;
    let mut seed = 0x5EED_1234_ABCD_0001u64;
    let mut random_axes = Vec::with_capacity(n_total);
    for _ in 0..n_total {
        let axis = random_unit_axis(&mut seed);
        random_axes.push(RandomAxisValue {
            axis,
            m: bishop_hill_taylor_factor_axis(axis, 1.0, slips).taylor_factor,
        });
    }
    let random_average = random_axes.iter().map(|a| a.m).sum::<f64>() / n_total as f64;
    random_axes.truncate(64);
    TaylorFile {
        label: label.to_string(),
        slip_family,
        n_random_total: n_total,
        random_average,
        axes,
        random_axes,
    }
}
// ---------------------------------------------------------------------
// 2. FCC single-crystal uniaxial tension (Phase 2 golden data).
// ---------------------------------------------------------------------

#[derive(Serialize)]
struct CpTensionStep {
    eps_11: f64,
    sigma_11: f64,
    sum_abs_slip_rates: f64,
    sum_accumulated_shear: f64,
    mean_crss: f64,
}

#[derive(Serialize)]
struct CpTensionFile {
    steps: Vec<CpTensionStep>,
}

fn cp_tension_dataset() -> CpTensionFile {
    let elastic = SymmetricFourthOrder::cubic(168.4e9, 121.4e9, 75.4e9);
    let model = CrystalPlasticityModel::from_crystal_structure(
        CrystalStructure::FCC,
        Hardening::Voce(VoceHardening::uniform(VoceParams {
            tau_0: 30.0e6,
            tau_s: 60.0e6,
            theta_0: 500.0e6,
            gamma_c: 0.05,
        })),
        Default::default(),
        elastic,
    )
    .unwrap();
    let mut state = HardeningState::from_hardening(&model.slip_systems, &model.hardening_law);
    let mut steps = Vec::with_capacity(8);
    for k in 1..=8 {
        let eps_inc = tpt_math_linalg_fixed::Vec6::new(0.001 * k as f64, 0.0, 0.0, 0.0, 0.0, 0.0);
        let update = solve_increment_single_point(&model, &mut state, &eps_inc).unwrap();
        steps.push(CpTensionStep {
            eps_11: 0.001 * k as f64,
            sigma_11: update.stress.data[0],
            sum_abs_slip_rates: update.slip_rates.iter().map(|g| g.abs()).sum(),
            sum_accumulated_shear: state.accumulated_shear.iter().sum(),
            mean_crss: state.crss.iter().sum::<f64>() / state.crss.len() as f64,
        });
    }
    CpTensionFile { steps }
}

// ---------------------------------------------------------------------
// 3. Phase-field goldens (Phase 3 golden data).
// ---------------------------------------------------------------------

#[derive(Serialize)]
struct SpinodalFile {
    grid_n: usize,
    n_steps: usize,
    mean_initial: f64,
    energy_initial: f64,
    energy_final: f64,
    interface_area_initial: f64,
    interface_area_final: f64,
    mean_final: f64,
}

fn spinodal_dataset() -> SpinodalFile {
    let n = 64usize;
    let grid = Grid2D::new(n, n, 1.0);
    let mut c = vec![0.5_f64; n * n];
    for (idx, v) in c.iter_mut().enumerate() {
        let i = idx / n;
        let j = idx % n;
        let phi = 2.0 * std::f64::consts::PI * ((i as f64) * 0.05 + (j as f64) * 0.07);
        *v = 0.5 + 0.05 * phi.sin() + 0.03 * (2.0 * phi).cos();
    }
    let mean_initial = c.iter().sum::<f64>() / c.len() as f64;
    let mut solver = PhaseFieldSolver::new_cahn_hilliard(
        grid,
        0.001,
        1.0,
        1.0,
        BulkEnergy::RegularSolution(RegularSolutionParams {
            omega: 4.0,
            rt: 1.0,
        }),
        &c,
    )
    .unwrap();
    let snap0 = solver.snapshot();
    let n_steps = 300usize;
    solver.step_many(n_steps).unwrap();
    let snap = solver.snapshot();
    let mean_final = solver.concentration.iter().sum::<f64>() / solver.concentration.len() as f64;
    SpinodalFile {
        grid_n: n,
        n_steps,
        mean_initial,
        energy_initial: snap0.free_energy,
        energy_final: snap.free_energy,
        interface_area_initial: snap0.interface_area,
        interface_area_final: snap.interface_area,
        mean_final,
    }
}
#[derive(Serialize)]
struct DendriteFile {
    grid_n: usize,
    n_steps: usize,
    solid_fraction_initial: f64,
    solid_fraction_final: f64,
    tip_velocity: f64,
    primary_arm_spacing: f64,
    secondary_arm_spacing: f64,
}

fn dendrite_dataset() -> DendriteFile {
    let n = 64usize;
    let grid = Grid2D::new(n, n, 1.0);
    let mut eta = vec![0.0_f64; n * n];
    let cx = n as f64 / 2.0;
    let cy = n as f64 / 2.0;
    let seed_radius = 4.0_f64;
    for i in 0..n {
        for j in 0..n {
            let d = ((i as f64 - cy).powi(2) + (j as f64 - cx).powi(2)).sqrt();
            eta[i * n + j] = if d < seed_radius { 1.0 } else { 0.0 };
        }
    }
    let anisotropy = AnisotropyModel {
        strength: 0.04,
        mode: AnisotropyMode::Cubic4Fold,
    };
    let mut solver = SolidificationSolver::new(
        grid,
        &eta,
        anisotropy,
        BulkEnergy::DoubleWell { well_depth: 1.0 },
        0.005,
        1.0,
        1.0,
    )
    .unwrap();
    let snap0 = solver.snapshot();
    let n_steps = 100usize;
    for _ in 0..n_steps {
        solver.step().unwrap();
    }
    let snap = solver.snapshot();
    DendriteFile {
        grid_n: n,
        n_steps,
        solid_fraction_initial: snap0.solid_fraction,
        solid_fraction_final: snap.solid_fraction,
        tip_velocity: snap.tip_velocity,
        primary_arm_spacing: snap.primary_arm_spacing,
        secondary_arm_spacing: snap.secondary_arm_spacing,
    }
}

// ---------------------------------------------------------------------
// 4. Degradation goldens (GTN void growth + fatigue crack initiation).
// ---------------------------------------------------------------------

#[derive(Serialize)]
struct YieldPoint {
    sigma_eq_over_y: f64,
    sigma_h_over_y: f64,
    f: f64,
    value: f64,
}

#[derive(Serialize)]
struct GtnFile {
    params: &'static str,
    f_initial: f64,
    f_trajectory: Vec<f64>,
    yield_grid: Vec<YieldPoint>,
}

fn gtn_dataset() -> GtnFile {
    let params = GtnParams::classic();
    let mut f = 0.001;
    let mut f_trajectory = Vec::with_capacity(12);
    for _ in 0..12 {
        f = update_porosity(f, 0.01, 0.05, &params);
        f_trajectory.push(f);
    }
    let mut yield_grid = Vec::new();
    for seq in [0.5, 0.75, 1.0, 1.25, 1.5] {
        for sh in [-1.0, -0.5, 0.0, 0.5, 1.0] {
            for fv in [0.001, 0.05, 0.2] {
                yield_grid.push(YieldPoint {
                    sigma_eq_over_y: seq,
                    sigma_h_over_y: sh,
                    f: fv,
                    value: gtn_yield_function(seq, sh, 1.0, fv, &params),
                });
            }
        }
    }
    GtnFile {
        params: "classic",
        f_initial: 0.001,
        f_trajectory,
        yield_grid,
    }
}
#[derive(Serialize)]
struct FatigueCase {
    criterion: String,
    cycles_simulated: f64,
    delta_gamma_per_cycle: f64,
    accumulated_shear: Vec<Vec<f64>>,
    cycles_to_initiation: f64,
    critical_grain: usize,
    critical_fip: f64,
}

#[derive(Serialize)]
struct FatigueFile {
    cases: Vec<FatigueCase>,
}

fn fatigue_dataset() -> FatigueFile {
    let mut seed = 0xC0FFEE_2026u64;
    let mut shear_history: Vec<Vec<f64>> = (0..6)
        .map(|_| (0..12).map(|_| 0.001 + 0.009 * lcg(&mut seed)).collect())
        .collect();
    shear_history[3] = vec![0.04; 12];
    shear_history[4] = vec![0.0; 12];

    let criteria: [(&str, FatigueCriterion); 3] = [
        ("findley k=0.3", FatigueCriterion::Findley { k: 0.3 }),
        (
            "fatemi-socie",
            FatigueCriterion::FatemiSocie {
                sigma_y: 200.0,
                k: 0.8,
            },
        ),
        (
            "crystallographic gamma_c=0.5",
            FatigueCriterion::CrystallographicSlip {
                critical_accumulated_shear: 0.5,
            },
        ),
    ];
    let cases = criteria
        .iter()
        .map(|(name, criterion)| {
            let r = predict_crack_initiation(&shear_history, criterion.clone(), 1000.0, 0.25);
            FatigueCase {
                criterion: name.to_string(),
                cycles_simulated: 1000.0,
                delta_gamma_per_cycle: 0.25,
                accumulated_shear: shear_history.clone(),
                cycles_to_initiation: r.cycles_to_initiation,
                critical_grain: r.critical_grain,
                critical_fip: r.critical_fip,
            }
        })
        .collect();
    FatigueFile { cases }
}

// ---------------------------------------------------------------------
// 5. Energy-materials golden (battery capacity fade).
// ---------------------------------------------------------------------

#[derive(Serialize)]
struct BatteryMilestone {
    cycle: usize,
    capacity_retention: f64,
    resistance_growth: f64,
}

#[derive(Serialize)]
struct BatteryFile {
    chemistry: &'static str,
    c_rate: f64,
    temperature_k: f64,
    num_cycles: usize,
    final_capacity_retention: f64,
    final_resistance_growth: f64,
    milestones: Vec<BatteryMilestone>,
}

fn battery_dataset() -> BatteryFile {
    let material = default_active_material(BatteryChemistry::Nmc811);
    let num_cycles = 1000usize;
    let curve = capacity_fade_curve(&material, 1.0, num_cycles, 298.0);
    let milestones = [1usize, 10, 100, 500, 1000]
        .into_iter()
        .map(|cycle| BatteryMilestone {
            cycle,
            capacity_retention: curve.capacity_retention[cycle - 1],
            resistance_growth: curve.resistance_growth[cycle - 1],
        })
        .collect();
    BatteryFile {
        chemistry: "NMC811",
        c_rate: 1.0,
        temperature_k: 298.0,
        num_cycles,
        final_capacity_retention: *curve.capacity_retention.last().unwrap(),
        final_resistance_growth: *curve.resistance_growth.last().unwrap(),
        milestones,
    }
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct ElectrodeGolden {
    equilibrium_potential: f64,
    exchange_current_density: f64,
    alpha_a: f64,
    alpha_c: f64,
    n: f64,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct CorrosionRateGolden {
    current_density: f64,
    corrosion_potential: f64,
    penetration_rate_mm_per_yr: f64,
    mass_loss_rate_g_per_m2_day: f64,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct PolarizationGolden {
    e_min: f64,
    e_max: f64,
    num_points: usize,
    potentials: Vec<f64>,
    currents: Vec<f64>,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
struct CorrosionFile {
    system: String,
    ph: f64,
    temperature_k: f64,
    anode: ElectrodeGolden,
    cathode: ElectrodeGolden,
    corrosion_fe: CorrosionRateGolden,
    corrosion_ti: CorrosionRateGolden,
    polarization: PolarizationGolden,
}

fn electrode_golden(e_eq: f64, i_0: f64, alpha_a: f64, alpha_c: f64, n: f64) -> ElectrodeGolden {
    ElectrodeGolden {
        equilibrium_potential: e_eq,
        exchange_current_density: i_0,
        alpha_a,
        alpha_c,
        n,
    }
}

fn corrosion_dataset() -> CorrosionFile {
    // Fe → Fe²⁺ + 2e⁻  vs  H₂ evolution, pH 0, 298 K (example default).
    let fe_anode = ElectrodeKinetics::from_alphas(-0.44, 1.0e-3, 0.5, 0.5, 2.0, 298.0);
    let h2_cathode = ElectrodeKinetics::from_alphas(0.0, 1.0e-1, 0.5, 0.5, 2.0, 298.0);
    let fe_model = CorrosionModel::new(fe_anode, h2_cathode, 0.0, 298.0);
    let cr_fe = corrosion_rate(&fe_model, 0.055_845, 2.0, 7874.0);

    let ti_anode = ElectrodeKinetics::from_alphas(-0.86, 1.0e-7, 0.5, 0.5, 3.0, 298.0);
    let ti_model = CorrosionModel::new(ti_anode, h2_cathode, 0.0, 298.0);
    let cr_ti = corrosion_rate(&ti_model, 0.047_867, 3.0, 4506.0);

    let pc = polarization_curve(&fe_model, (-0.6, 0.2), 33, PolarizationBranch::Net);

    CorrosionFile {
        system: "Fe / H2 (pH 0, 298 K)".to_string(),
        ph: 0.0,
        temperature_k: 298.0,
        anode: electrode_golden(-0.44, 1.0e-3, 0.5, 0.5, 2.0),
        cathode: electrode_golden(0.0, 1.0e-1, 0.5, 0.5, 2.0),
        corrosion_fe: CorrosionRateGolden {
            current_density: cr_fe.current_density,
            corrosion_potential: cr_fe.corrosion_potential,
            penetration_rate_mm_per_yr: cr_fe.penetration_rate_mm_per_yr,
            mass_loss_rate_g_per_m2_day: cr_fe.mass_loss_rate_g_per_m2_day,
        },
        corrosion_ti: CorrosionRateGolden {
            current_density: cr_ti.current_density,
            corrosion_potential: cr_ti.corrosion_potential,
            penetration_rate_mm_per_yr: cr_ti.penetration_rate_mm_per_yr,
            mass_loss_rate_g_per_m2_day: cr_ti.mass_loss_rate_g_per_m2_day,
        },
        polarization: PolarizationGolden {
            e_min: -0.6,
            e_max: 0.2,
            num_points: pc.potentials.len(),
            potentials: pc.potentials,
            currents: pc.currents,
        },
    }
}

// ---------------------------------------------------------------------

fn main() {
    write_json(
        "taylor-factor/fcc.json",
        &taylor_dataset("fcc", "{111}<110>", None),
    );
    let bcc = CrystalStructure::BCC.slip_systems();
    write_json(
        "taylor-factor/bcc.json",
        &taylor_dataset("bcc", "{110}<111>", Some(&bcc)),
    );
    write_json(
        "crystal-plasticity/fcc-single-crystal-tension.json",
        &cp_tension_dataset(),
    );
    write_json(
        "phase-field/spinodal-decomposition.json",
        &spinodal_dataset(),
    );
    write_json(
        "phase-field/dendritic-solidification.json",
        &dendrite_dataset(),
    );
    write_json("degradation/gtn-void-growth.json", &gtn_dataset());
    write_json(
        "degradation/fatigue-crack-initiation.json",
        &fatigue_dataset(),
    );
    write_json(
        "energy-materials/battery-degradation.json",
        &battery_dataset(),
    );
    write_json(
        "energy-materials/corrosion-polarization.json",
        &corrosion_dataset(),
    );
    println!("golden datasets regenerated");
}
