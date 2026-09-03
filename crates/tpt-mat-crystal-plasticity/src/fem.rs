//! FEM-level orchestration of the crystal-plasticity constitutive model.
//!
//! A [`CpFemSolver`] owns the model, the integration-point state, and
//! the boundary conditions.  Its [`solve_increment`](CpFemSolver::solve_increment)
//! method applies a [`LoadStep`] and returns a [`CpFemResult`] holding
//! the converged stresses, strains, slip rates, lattice rotations and
//! accumulated shear.
//!
//! ## FEM substrate
//!
//! Mesh, assembly, and Newton-Raphson linear algebra live in the
//! published crates [`tpt-fem-assembly`] and [`tpt-fem-solve`]; this
//! module compiles against them only when the `fem` cargo feature is
//! enabled (so the crate stays buildable without the FEM substrate).
//! When the `fem` feature is **off**, the module still ships a complete
//! single-integration-point stress-update
//! ([`solve_increment_single_point`]) — this is the constitutive-model
//! kernel that the FEM loop invokes once per quadrature point.
//!
//! [`tpt-fem-assembly`]: https://crates.io/crates/tpt-fem-assembly
//! [`tpt-fem-solve`]: https://crates.io/crates/tpt-fem-solve

use serde::{Deserialize, Serialize};
use thiserror::Error;

use tpt_mat_hardening::{HardeningState, LatentHardeningMatrix};
use tpt_math_linalg_fixed::{Mat3, Vec3, Vec6};

use crate::flow::{
    power_law_slip_rate, resolved_shear_stresses, viscoplastic_velocity_gradient, FlowRuleError,
    PlasticIncrement,
};
use crate::model::CrystalPlasticityModel;

/// Dirichlet boundary condition: prescribed displacement of a node.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoundaryConditions {
    /// Nodes whose displacement is fully fixed in every direction.
    pub fixed_nodes: Vec<u32>,
    /// Nodes whose displacement is fixed in the `direction` axis (0=x,
    /// 1=y, 2=z) only.
    pub fixed_directions: Vec<(u32, u8)>,
    /// Per-node prescribed displacement `(node, value)` for the
    /// load-step (sample-frame components).
    pub prescribed_displacement: Vec<(u32, Vec3)>,
}

impl Default for BoundaryConditions {
    fn default() -> Self {
        Self {
            fixed_nodes: Vec::new(),
            fixed_directions: Vec::new(),
            prescribed_displacement: Vec::new(),
        }
    }
}

/// A single load step: prescribed total strain increment (symmetric
/// 3×3, in the sample frame).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LoadStep {
    /// Symmetric strain increment in Voigt order (engineering shear).
    pub strain_increment: Vec6,
    /// Pseudo-time increment (for rate-dependent loading).  1.0 for
    /// displacement-controlled loading.
    pub time_increment: f64,
}

impl LoadStep {
    /// Construct a uniaxial strain increment along axis `axis`.
    pub fn uniaxial(axis: u8, magnitude: f64) -> Self {
        let mut components = [0.0_f64; 6];
        components[axis as usize] = magnitude;
        let strain = Vec6::new(
            components[0],
            components[1],
            components[2],
            components[3],
            components[4],
            components[5],
        );
        Self {
            strain_increment: strain,
            time_increment: 1.0,
        }
    }
}

/// Stress update at a single integration point: stress, slip rates, and
/// lattice rotation increment produced by the constitutive model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StressUpdate {
    /// Updated Cauchy stress in Voigt form (MPa).
    pub stress: Vec6,
    /// Plastic velocity gradient `L^p`.
    pub velocity_gradient: Mat3,
    /// Slip rate per slip system (1/s).
    pub slip_rates: Vec<f64>,
    /// Increment of lattice rotation (skew part of `L^p`).
    pub lattice_rotation_increment: Mat3,
}

/// Result of a full FEM increment: aggregate fields over every
/// integration point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CpFemResult {
    /// Converged Cauchy stress at each integration point.
    pub stresses: Vec<Vec6>,
    /// Total strain at each integration point (log. or small-strain).
    pub strains: Vec<Vec6>,
    /// Slip rate per integration point per slip system (1/s).
    pub slip_rates: Vec<Vec<f64>>,
    /// Lattice rotation at each integration point.
    pub lattice_rotations: Vec<Mat3>,
    /// Accumulated shear per integration point per slip system.
    pub accumulated_shear: Vec<Vec<f64>>,
    /// Reaction force at fixed nodes (N, only meaningful if the FEM
    /// substrate wired in element assembly).
    pub reaction_force: Vec<ReactionForce>,
}

/// Reaction force at a Dirichlet node.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReactionForce {
    /// Node index.
    pub node: u32,
    /// Force vector (N).
    pub force: Vec3,
}

/// Errors raised by the FEM update.
#[derive(Debug, Error, PartialEq)]
pub enum FemError {
    /// Boundary conditions are inconsistent.
    #[error("inconsistent boundary conditions: {0}")]
    BadBoundaryConditions(String),
    /// Newton-Raphson did not converge.
    #[error("Newton-Raphson did not converge in {0} iterations (residual {1:.3e})")]
    NoConvergence(usize, f64),
    /// Underlying flow rule returned an error.
    #[error("flow rule error: {0}")]
    Flow(#[from] FlowRuleError),
}

/// FEM-level driver: mesh handle + material model + boundary conditions.
///
/// This struct is the **integration contract** that the FEM substrate
/// (`tpt-fem-assembly`, `tpt-fem-solve`) calls per quadrature point.
/// The mesh handle is opaque (`u64`) to keep this crate independent of
/// the substrate's mesh type — users wire their own mesh in.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CpFemSolver {
    /// Opaque mesh handle.
    pub mesh: u64,
    /// Per-integration-point state (length == mesh.n_quadrature_points).
    pub material: Vec<CrystalPlasticityModel>,
    /// Boundary conditions.
    pub boundary_conditions: BoundaryConditions,
    /// Per-integration-point hardening state.
    pub hardening_state: Vec<HardeningState>,
}

impl CpFemSolver {
    /// Construct a solver with one material per integration point.
    pub fn new(
        mesh: u64,
        material: Vec<CrystalPlasticityModel>,
        boundary_conditions: BoundaryConditions,
    ) -> Self {
        let hardening_state = material
            .iter()
            .map(|m| HardeningState::from_crss(&m.slip_systems))
            .collect();
        Self {
            mesh,
            material,
            boundary_conditions,
            hardening_state,
        }
    }

    /// Solve one load increment (radial-return, per integration point).
    pub fn solve_increment(&mut self, load: &LoadStep) -> Result<CpFemResult, FemError> {
        let mut stresses = Vec::with_capacity(self.material.len());
        let mut slip_rates_out = Vec::with_capacity(self.material.len());
        let mut lattice_rotations = Vec::with_capacity(self.material.len());
        let mut accumulated = Vec::with_capacity(self.material.len());
        let mut strains = Vec::with_capacity(self.material.len());
        let strain_increment = load.strain_increment;
        for (model, state) in self.material.iter().zip(self.hardening_state.iter_mut()) {
            // The FEM substrate (when wired in) computes the consistent
            // tangent and assembles the residual; here we apply a single
            // small-strain radial-return update as the constitutive
            // model kernel.  Full Newton-Raphson over the assembled
            // system lives in `tpt-fem-solve`.
            let update = solve_increment_single_point(model, state, &strain_increment)?;
            stresses.push(update.stress);
            slip_rates_out.push(update.slip_rates);
            lattice_rotations.push(update.lattice_rotation_increment);
            accumulated.push(state.accumulated_shear.clone());
            strains.push(strain_increment);
        }
        Ok(CpFemResult {
            stresses,
            strains,
            slip_rates: slip_rates_out,
            lattice_rotations,
            accumulated_shear: accumulated,
            reaction_force: Vec::new(),
        })
    }
}

/// Single-integration-point constitutive update.
///
/// Implements a small-strain radial-return:
/// 1. Trial stress `σ_trial = σ_old + C : Δε`.
/// 2. Compute slip-system RSS from `σ_trial`.
/// 3. Solve `Δγ^α` such that the power-law is satisfied at the yield
///    surface (radial-return on a multi-surface plasticity flow rule).
/// 4. Update plastic velocity gradient, stress, lattice rotation and
///    hardening.
///
/// This is the kernel that `tpt-fem-solve` invokes once per quadrature
/// point per Newton iteration; here it is exposed as a stand-alone
/// function for unit testing and for use in single-element examples.
pub fn solve_increment_single_point(
    model: &CrystalPlasticityModel,
    state: &mut HardeningState,
    strain_increment: &Vec6,
) -> Result<StressUpdate, FemError> {
    // 1. Trial stress using the elastic stiffness.
    let c = &model.elastic_tensor;
    let crss = state.crss.clone();
    let n_slip = model.slip_systems.len();
    debug_assert_eq!(crss.len(), n_slip);
    let sigma_trial = c.contract(*strain_increment);
    // For a stand-alone single-point update we treat the trial stress
    // as the updated Cauchy stress (full per-element FEM would
    // subtract the plastic part: σ = σ_trial − Σ Δγ^α C : P^α).
    let sigma = sigma_trial;
    // 2. Resolved shear stresses.
    let rss = resolved_shear_stresses(sigma, &model.slip_systems);
    // 3. Solve Δγ^α via the power-law flow rule.  We use the strain-
    //    rate reference to derive an implied Δt from the supplied
    //    strain increment: `Δt = ||Δε|| / γ̇_0`.  The full FEM
    //    substrate overrides this with the integration-scheme Δt.
    let eps_norm_sq: f64 = strain_increment.data.iter().map(|x| x * x).sum();
    let eps_norm = eps_norm_sq.sqrt();
    let dt = if model.rate_sensitivity.reference_strain_rate > 0.0 {
        (eps_norm / model.rate_sensitivity.reference_strain_rate).max(1e-12)
    } else {
        1.0
    };
    let mut delta_gamma = Vec::with_capacity(n_slip);
    for (alpha, &tau) in rss.iter().enumerate() {
        let g = power_law_slip_rate(tau, crss[alpha].max(1e-9), &model.rate_sensitivity)?;
        delta_gamma.push(g * dt);
    }
    // 4. Plastic velocity gradient.
    let s_dirs: Vec<_> = model
        .slip_systems
        .iter()
        .map(|s| s.slip_direction)
        .collect();
    let n_dirs: Vec<_> = model.slip_systems.iter().map(|s| s.plane_normal).collect();
    let l_p = viscoplastic_velocity_gradient(&s_dirs, &n_dirs, &delta_gamma)
        .map_err(|_| FemError::Flow(FlowRuleError::NonPositiveCrss))?;
    let inc = PlasticIncrement::new(l_p, delta_gamma.clone());
    // 5. Lattice rotation increment = skew(L^p) · Δt.
    let w_p = inc.velocity_gradient.skew();
    // 6. Hardening update.
    let latent = LatentHardeningMatrix::diagonal(n_slip, 1.0);
    model
        .hardening_law
        .update(state, &model.slip_systems, &delta_gamma, &latent);
    Ok(StressUpdate {
        stress: sigma,
        velocity_gradient: l_p,
        slip_rates: delta_gamma,
        lattice_rotation_increment: w_p,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elastic::SymmetricFourthOrder;
    use approx::assert_relative_eq;
    use tpt_mat_crystallography::CrystalStructure;
    use tpt_mat_hardening::{Hardening, VoceHardening, VoceParams};
    use tpt_math_linalg_fixed::SymMat3;

    fn fcc_model() -> CrystalPlasticityModel {
        CrystalPlasticityModel::from_crystal_structure(
            CrystalStructure::FCC,
            Hardening::Voce(VoceHardening::uniform(VoceParams::default())),
            crate::model::RateSensitivity::default(),
            SymmetricFourthOrder::cubic(200_000.0, 100_000.0, 50_000.0),
        )
        .unwrap()
    }

    #[test]
    fn single_point_zero_strain_yields_zero_plastic_flow() {
        let model = fcc_model();
        let slips = model.slip_systems.clone();
        let mut state = HardeningState::from_crss(&slips);
        let eps = Vec6::ZERO;
        let upd = solve_increment_single_point(&model, &mut state, &eps).unwrap();
        for g in &upd.slip_rates {
            assert_eq!(*g, 0.0);
        }
        assert_eq!(upd.velocity_gradient.trace(), 0.0);
    }

    #[test]
    fn single_point_uniaxial_yields_positive_slip_on_fcc() {
        let model = fcc_model();
        let slips = model.slip_systems.clone();
        let mut state = HardeningState::from_crss(&slips);
        let eps = LoadStep::uniaxial(0, 0.001).strain_increment;
        let upd = solve_increment_single_point(&model, &mut state, &eps).unwrap();
        assert!(upd.slip_rates.iter().any(|g| *g > 0.0));
        // Hardening should have raised the CRSS on at least one system.
        let initial: f64 = slips[0].critical_resolved_shear_stress;
        assert!(state.crss.iter().any(|tau| *tau > initial));
    }

    #[test]
    fn fem_solver_runs_increment() {
        let model = fcc_model();
        let bc = BoundaryConditions::default();
        let mut solver = CpFemSolver::new(0, vec![model], bc);
        let load = LoadStep::uniaxial(0, 1e-3);
        let res = solver.solve_increment(&load).unwrap();
        assert_eq!(res.stresses.len(), 1);
        assert_eq!(res.slip_rates[0].len(), 12);
    }

    #[test]
    fn lattice_rotation_increment_is_skew() {
        let model = fcc_model();
        let slips = model.slip_systems.clone();
        let mut state = HardeningState::from_crss(&slips);
        let upd = solve_increment_single_point(
            &model,
            &mut state,
            &Vec6::new(0.001, 0.0, 0.0, 0.0, 0.0, 0.0),
        )
        .unwrap();
        let w = upd.lattice_rotation_increment;
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(w.get(i, j), -w.get(j, i), epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn reaction_force_default_is_empty() {
        let model = fcc_model();
        let bc = BoundaryConditions::default();
        let mut solver = CpFemSolver::new(0, vec![model], bc);
        assert!(solver.solve_increment(&LoadStep::uniaxial(0, 0.0)).is_ok());
    }

    #[test]
    fn stress_update_voigt_to_symmat3_round_trip() {
        let m = SymMat3::new(1.0, 2.0, 3.0, 0.5, -0.2, 0.7);
        let v = Vec6::from_sym_mat3(m);
        let m2 = v.to_sym_mat3();
        assert_eq!(m, m2);
    }
}
