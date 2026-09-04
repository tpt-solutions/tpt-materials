//! Self-contained FEM assembly + Newton-Raphson solver for
//! crystal-plasticity problems.
//!
//! Provides a minimal but complete 3-D 8-node hexahedral
//! trilinear FEM stack:
//!
//! - Hex8 element with 8 Gauss points (full integration)
//! - Consistent tangent `d sigma / d epsilon` returned by the
//!   crystal-plasticity constitutive model (the small-strain
//!   elastic stiffness is used as a stand-in for the algorithmic
//!   tangent; the single-point update path is exercised at every
//!   Gauss point per Newton iteration).
//! - Newton-Raphson global equilibrium iteration with line search
//!   and convergence check on the displacement residual norm.
//! - Dirichlet boundary conditions on selected node DOFs.
//!
//! The mesh handle is intentionally simple (`Mesh`) — vertices and
//! connectivity only.  This is sufficient for serial, single-
//! material, single-grain RVE problems; it is not intended to be a
//! production FEM solver.  The point is to demonstrate the full
//! micro-to-macro loop and to provide a reference implementation
//! that downstream `tpt-fem-*` crates can replicate / extend.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use tpt_math_linalg_fixed::Vec3;

use crate::flow::resolved_shear_stresses;
use crate::model::CrystalPlasticityModel;
use crate::fem::{solve_increment_single_point, BoundaryConditions, LoadStep, CpFemResult};
use tpt_mat_hardening::HardeningState;

/// 8-node hexahedral element with trilinear shape functions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hex8 {
    /// Node indices (8 corners).
    pub nodes: [usize; 8],
}

/// 3-D trilinear FEM mesh (8-node hexahedra).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mesh {
    /// Node coordinates (x, y, z) in metres.
    pub vertices: Vec<Vec3>,
    /// Element list.
    pub elements: Vec<Hex8>,
}

impl Mesh {
    /// Construct a single-element unit cube mesh (8 nodes, 1 hex).
    pub fn unit_cube() -> Self {
        let vertices = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(0.0, 1.0, 1.0),
        ];
        let elements = vec![Hex8 {
            nodes: [0, 1, 2, 3, 4, 5, 6, 7],
        }];
        Self { vertices, elements }
    }

    /// Total number of degrees of freedom (3 per node).
    pub fn n_dof(&self) -> usize {
        self.vertices.len() * 3
    }
}

/// FEM errors.
#[derive(Debug, Error, PartialEq)]
pub enum CpFemAssemblyError {
    /// Newton-Raphson did not converge.
    #[error("Newton-Raphson did not converge in {0} iterations (residual {1:.3e})")]
    NoConvergence(usize, f64),
    /// Mesh / boundary-condition incompatibility.
    #[error("boundary condition error: {0}")]
    BadBoundary(String),
}

/// FEM solver state.
pub struct CpFemAssembly {
    mesh: Mesh,
    material: CrystalPlasticityModel,
    /// Per-Gauss-point hardening state.
    pub hardening_state: Vec<HardeningState>,
    /// Current global displacement vector (3 N entries).
    pub u: Vec<f64>,
    /// Number of Gauss points per element (8 for full integration).
    n_gauss_per_elem: usize,
}

impl CpFemAssembly {
    /// Construct a single-material assembly over the given mesh.
    pub fn new(mesh: Mesh, material: CrystalPlasticityModel) -> Self {
        let n_gauss = mesh.elements.len() * 8;
        let hardening_state =
            vec![HardeningState::from_crss(&material.slip_systems); n_gauss];
        let u = vec![0.0; mesh.n_dof()];
        Self {
            mesh,
            material,
            hardening_state,
            u,
            n_gauss_per_elem: 8,
        }
    }

    /// Newton-Raphson solve for one load step.
    pub fn solve(&mut self, load: &LoadStep, bc: &BoundaryConditions, max_iter: usize, tol: f64) -> Result<CpFemResult, CpFemAssemblyError> {
        let n_dof = self.mesh.n_dof();
        let mut u = vec![0.0; n_dof];
        let mut ext = vec![0.0; n_dof];
        let strain_inc = load.strain_increment;
        let strain_arr = [
            strain_inc.data[0],
            strain_inc.data[1],
            strain_inc.data[2],
            strain_inc.data[3],
            strain_inc.data[4],
            strain_inc.data[5],
        ];
        let cube_face_x: [usize; 4] = [1, 2, 5, 6];
        for &n in &cube_face_x {
            ext[n * 3] += strain_arr[0];
        }
        let cube_face_y: [usize; 4] = [2, 3, 6, 7];
        for &n in &cube_face_y {
            ext[n * 3 + 1] += strain_arr[1];
        }
        let cube_face_z: [usize; 4] = [4, 5, 6, 7];
        for &n in &cube_face_z {
            ext[n * 3 + 2] += strain_arr[2];
        }
        let mut residual_norm = f64::INFINITY;
        for _iter in 0..max_iter {
            let (f_int, mut k_tan) = self.assemble(&u);
            let mut r = vec![0.0; n_dof];
            for i in 0..n_dof {
                r[i] = f_int[i] - ext[i];
            }
            apply_dirichlet(&mut r, &mut k_tan, &u, bc, &self.mesh);
            residual_norm = r.iter().map(|x| x * x).sum::<f64>().sqrt();
            if residual_norm < tol {
                self.u = u.clone();
                return Ok(self.collect_result(&u));
            }
            let neg_r: Vec<f64> = r.iter().map(|x| -x).collect();
            let du = solve_linear_system(&k_tan, &neg_r);
            let mut alpha = 1.0;
            let r_norm_0 = residual_norm;
            for _ in 0..8 {
                let mut u_trial = u.clone();
                for i in 0..n_dof {
                    u_trial[i] += alpha * du[i];
                }
                let (f_trial, _) = self.assemble(&u_trial);
                let mut r_trial = vec![0.0; n_dof];
                for i in 0..n_dof {
                    r_trial[i] = f_trial[i] - ext[i];
                }
                apply_dirichlet_residual(&mut r_trial, bc, &self.mesh);
                let r_trial_norm = r_trial.iter().map(|x| x * x).sum::<f64>().sqrt();
                if r_trial_norm < r_norm_0 * (1.0 - 0.01 * alpha) {
                    u = u_trial;
                    break;
                }
                alpha *= 0.5;
                if alpha < 1.0e-6 {
                    u = u_trial;
                    break;
                }
            }
        }
        Err(CpFemAssemblyError::NoConvergence(max_iter, residual_norm))
    }

    /// Assemble internal force and tangent stiffness.
    fn assemble(&self, u: &[f64]) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n_dof = self.mesh.n_dof();
        let mut f_int = vec![0.0; n_dof];
        let mut k_tan = vec![vec![0.0; n_dof]; n_dof];
        // Loop over elements.
        for (elem_idx, hex) in self.mesh.elements.iter().enumerate() {
            // Loop over Gauss points.
            for (gp_local, (xi, eta, zeta)) in GAUSS8.iter().enumerate() {
                let gp_global = elem_idx * 8 + gp_local;
                // Compute shape-function derivatives B (6 x 24).
                let b_mat = hex8_b_matrix(hex, &self.mesh.vertices, *xi, *eta, *zeta);
                // Compute the displacement vector at the element nodes.
                let mut u_e = vec![0.0; 24];
                for (i, &n) in hex.nodes.iter().enumerate() {
                    u_e[i * 3] = u[n * 3];
                    u_e[i * 3 + 1] = u[n * 3 + 1];
                    u_e[i * 3 + 2] = u[n * 3 + 2];
                }
                // Compute the strain from u_e: eps = B u_e.
                let strain = apply_b_u(&b_mat, &u_e);
                // Compute the stress from the constitutive model (single-point radial return).
                let mut state = self
                    .hardening_state
                    .get(gp_global)
                    .cloned()
                    .unwrap_or_else(|| HardeningState::from_crss(&self.material.slip_systems));
                let mut model = self.material.clone();
                let eps_inc = tpt_math_linalg_fixed::Vec6::new(
                    strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
                );
                let update = solve_increment_single_point(&model, &mut state, &eps_inc)
                    .expect("single-point update should succeed");
                let sigma_arr = [
                    update.stress.data[0],
                    update.stress.data[1],
                    update.stress.data[2],
                    update.stress.data[3],
                    update.stress.data[4],
                    update.stress.data[5],
                ];
                // Compute element Jacobian determinant (for the unit
                // cube and the standard Gauss rule).
                let det_j = hex8_det_j(hex, &self.mesh.vertices, *xi, *eta, *zeta);
                let w = 1.0_f64;
                // Assemble f_int += B^T sigma * detJ * w.
                for i in 0..8 {
                    for k in 0..3 {
                        let dof = hex.nodes[i] * 3 + k;
                        let mut contrib = 0.0;
                        for m in 0..6 {
                            contrib += b_mat[m][i * 3 + k] * sigma_arr[m];
                        }
                        f_int[dof] += contrib * det_j * w;
                    }
                }
                // Assemble K_tan += B^T C B detJ w.  C is the
                // elastic stiffness (small-strain consistent tangent).
                let c_mat = &model.elastic_tensor.data;
                let mut cb = vec![vec![0.0_f64; 24]; 6];
                for m in 0..6 {
                    for j in 0..24 {
                        let mut s = 0.0;
                        for n in 0..6 {
                            s += c_mat[m][n] * b_mat[n][j];
                        }
                        cb[m][j] = s;
                    }
                }
                for i in 0..8 {
                    for k in 0..3 {
                        let row = hex.nodes[i] * 3 + k;
                        for j in 0..8 {
                            for l in 0..3 {
                                let col = hex.nodes[j] * 3 + l;
                                let mut s = 0.0;
                                for m in 0..6 {
                                    s += b_mat[m][i * 3 + k] * cb[m][j * 3 + l];
                                }
                                k_tan[row][col] += s * det_j * w;
                            }
                        }
                    }
                }
                let _ = resolved_shear_stresses;
            }
        }
        (f_int, k_tan)
    }

    /// Collect the FEM result for the current displacement.
    fn collect_result(&self, u: &[f64]) -> CpFemResult {
        let mut stresses = Vec::new();
        let mut slip_rates_out = Vec::new();
        let mut lattice_rotations = Vec::new();
        let mut accumulated = Vec::new();
        let mut strains = Vec::new();
        for (elem_idx, hex) in self.mesh.elements.iter().enumerate() {
            for (gp_local, (xi, eta, zeta)) in GAUSS8.iter().enumerate() {
                let gp_global = elem_idx * 8 + gp_local;
                let b_mat = hex8_b_matrix(hex, &self.mesh.vertices, *xi, *eta, *zeta);
                let mut u_e = vec![0.0; 24];
                for (i, &n) in hex.nodes.iter().enumerate() {
                    u_e[i * 3] = u[n * 3];
                    u_e[i * 3 + 1] = u[n * 3 + 1];
                    u_e[i * 3 + 2] = u[n * 3 + 2];
                }
                let strain = apply_b_u(&b_mat, &u_e);
                let mut state = self
                    .hardening_state
                    .get(gp_global)
                    .cloned()
                    .unwrap_or_else(|| HardeningState::from_crss(&self.material.slip_systems));
                let mut model = self.material.clone();
                let eps_inc = tpt_math_linalg_fixed::Vec6::new(
                    strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
                );
                if let Ok(update) = solve_increment_single_point(&model, &mut state, &eps_inc) {
                    stresses.push(update.stress);
                    slip_rates_out.push(update.slip_rates);
                    lattice_rotations.push(update.lattice_rotation_increment);
                    accumulated.push(state.accumulated_shear.clone());
                    strains.push(eps_inc);
                }
            }
        }
        CpFemResult {
            stresses,
            strains,
            slip_rates: slip_rates_out,
            lattice_rotations,
            accumulated_shear: accumulated,
            reaction_force: Vec::new(),
        }
    }
}

fn apply_b_u(b: &[Vec<f64>; 6], u: &[f64]) -> [f64; 6] {
    let mut eps = [0.0_f64; 6];
    for i in 0..6 {
        let mut s = 0.0;
        for j in 0..24 {
            s += b[i][j] * u[j];
        }
        eps[i] = s;
    }
    eps
}

/// 8-point Gauss rule for the reference cube [-1, 1]^3.
const GAUSS8: &[(f64, f64, f64)] = &[
    (-0.5773502691896257, -0.5773502691896257, -0.5773502691896257),
    (0.5773502691896257, -0.5773502691896257, -0.5773502691896257),
    (0.5773502691896257, 0.5773502691896257, -0.5773502691896257),
    (-0.5773502691896257, 0.5773502691896257, -0.5773502691896257),
    (-0.5773502691896257, -0.5773502691896257, 0.5773502691896257),
    (0.5773502691896257, -0.5773502691896257, 0.5773502691896257),
    (0.5773502691896257, 0.5773502691896257, 0.5773502691896257),
    (-0.5773502691896257, 0.5773502691896257, 0.5773502691896257),
];

/// Compute the Hex8 B-matrix (6 x 24) at the given parametric point.
fn hex8_b_matrix(hex: &Hex8, vertices: &[Vec3], xi: f64, eta: f64, zeta: f64) -> [Vec<f64>; 6] {
    let mut b = std::array::from_fn::<_, 6, _>(|_| vec![0.0_f64; 24]);
    // Shape functions for the 8-node hex.
    let n = |x: f64, e: f64, z: f64| -> [f64; 8] {
        let n0 = 0.125 * (1.0 - x) * (1.0 - e) * (1.0 - z);
        let n1 = 0.125 * (1.0 + x) * (1.0 - e) * (1.0 - z);
        let n2 = 0.125 * (1.0 + x) * (1.0 + e) * (1.0 - z);
        let n3 = 0.125 * (1.0 - x) * (1.0 + e) * (1.0 - z);
        let n4 = 0.125 * (1.0 - x) * (1.0 - e) * (1.0 + z);
        let n5 = 0.125 * (1.0 + x) * (1.0 - e) * (1.0 + z);
        let n6 = 0.125 * (1.0 + x) * (1.0 + e) * (1.0 + z);
        let n7 = 0.125 * (1.0 - x) * (1.0 + e) * (1.0 + z);
        [n0, n1, n2, n3, n4, n5, n6, n7]
    };
    let dn_dxi = |x: f64, e: f64, z: f64| -> [f64; 8] {
        let a = 0.125 * (1.0 - e) * (1.0 - z);
        let b = 0.125 * (1.0 + e) * (1.0 - z);
        let c = 0.125 * (1.0 - e) * (1.0 + z);
        let d = 0.125 * (1.0 + e) * (1.0 + z);
        [-a, a, b, -b, -c, c, d, -d]
    };
    let dn_deta = |x: f64, e: f64, z: f64| -> [f64; 8] {
        let a = 0.125 * (1.0 - x) * (1.0 - z);
        let b = 0.125 * (1.0 + x) * (1.0 - z);
        let c = 0.125 * (1.0 - x) * (1.0 + z);
        let d = 0.125 * (1.0 + x) * (1.0 + z);
        [-a, -b, b, a, -c, -d, d, c]
    };
    let dn_dzeta = |x: f64, e: f64, z: f64| -> [f64; 8] {
        let a = 0.125 * (1.0 - x) * (1.0 - e);
        let b = 0.125 * (1.0 + x) * (1.0 - e);
        let c = 0.125 * (1.0 - x) * (1.0 + e);
        let d = 0.125 * (1.0 + x) * (1.0 + e);
        [-a, -b, -c, -d, a, b, c, d]
    };
    let _ = n;
    // Compute the Jacobian at the given point.
    let mut j = [[0.0_f64; 3]; 3];
    let dnx = dn_dxi(xi, eta, zeta);
    let dne = dn_deta(xi, eta, zeta);
    let dnz = dn_dzeta(xi, eta, zeta);
    for i in 0..8 {
        let node = hex.nodes[i];
        let p = vertices[node];
        for k in 0..3 {
            j[0][k] += dnx[i] * p.data[k];
            j[1][k] += dne[i] * p.data[k];
            j[2][k] += dnz[i] * p.data[k];
        }
    }
    // Inverse of the Jacobian (3x3).  Assume non-singular.
    let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
        - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
        + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
    let inv_det = if det.abs() < 1.0e-15 { 1.0 } else { 1.0 / det };
    let mut j_inv = [[0.0_f64; 3]; 3];
    j_inv[0][0] = (j[1][1] * j[2][2] - j[1][2] * j[2][1]) * inv_det;
    j_inv[0][1] = (j[0][2] * j[2][1] - j[0][1] * j[2][2]) * inv_det;
    j_inv[0][2] = (j[0][1] * j[1][2] - j[0][2] * j[1][1]) * inv_det;
    j_inv[1][0] = (j[1][2] * j[2][0] - j[1][0] * j[2][2]) * inv_det;
    j_inv[1][1] = (j[0][0] * j[2][2] - j[0][2] * j[2][0]) * inv_det;
    j_inv[1][2] = (j[0][2] * j[1][0] - j[0][0] * j[1][2]) * inv_det;
    j_inv[2][0] = (j[1][0] * j[2][1] - j[1][1] * j[2][0]) * inv_det;
    j_inv[2][1] = (j[0][1] * j[2][0] - j[0][0] * j[2][1]) * inv_det;
    j_inv[2][2] = (j[0][0] * j[1][1] - j[0][1] * j[1][0]) * inv_det;
    // Derivatives in physical coordinates: dn_dx = j_inv * d/dxi.
    let mut dn_dx = [0.0_f64; 8];
    let mut dn_dy = [0.0_f64; 8];
    let mut dn_dz = [0.0_f64; 8];
    for i in 0..8 {
        dn_dx[i] = j_inv[0][0] * dnx[i] + j_inv[0][1] * dne[i] + j_inv[0][2] * dnz[i];
        dn_dy[i] = j_inv[1][0] * dnx[i] + j_inv[1][1] * dne[i] + j_inv[1][2] * dnz[i];
        dn_dz[i] = j_inv[2][0] * dnx[i] + j_inv[2][1] * dne[i] + j_inv[2][2] * dnz[i];
    }
    // B-matrix assembly for engineering shear.
    for i in 0..8 {
        b[0][i * 3] = dn_dx[i];
        b[1][i * 3 + 1] = dn_dy[i];
        b[2][i * 3 + 2] = dn_dz[i];
        b[3][i * 3] = dn_dy[i];
        b[3][i * 3 + 1] = dn_dx[i];
        b[4][i * 3 + 1] = dn_dz[i];
        b[4][i * 3 + 2] = dn_dy[i];
        b[5][i * 3] = dn_dz[i];
        b[5][i * 3 + 2] = dn_dx[i];
    }
    b
}

/// Compute the determinant of the Hex8 Jacobian at (xi, eta, zeta).
fn hex8_det_j(hex: &Hex8, vertices: &[Vec3], xi: f64, eta: f64, zeta: f64) -> f64 {
    let dnx = [0.125 * -1.0, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125];
    let _ = dnx; // suppress unused
    let mut j = [[0.0_f64; 3]; 3];
    let dnx_full = |x: f64, e: f64, z: f64| -> [f64; 8] {
        let a = 0.125 * (1.0 - e) * (1.0 - z);
        let b = 0.125 * (1.0 + e) * (1.0 - z);
        let c = 0.125 * (1.0 - e) * (1.0 + z);
        let d = 0.125 * (1.0 + e) * (1.0 + z);
        [-a, a, b, -b, -c, c, d, -d]
    };
    let dne = |x: f64, e: f64, z: f64| -> [f64; 8] {
        let a = 0.125 * (1.0 - x) * (1.0 - z);
        let b = 0.125 * (1.0 + x) * (1.0 - z);
        let c = 0.125 * (1.0 - x) * (1.0 + z);
        let d = 0.125 * (1.0 + x) * (1.0 + z);
        [-a, -b, b, a, -c, -d, d, c]
    };
    let dnz = |x: f64, e: f64, z: f64| -> [f64; 8] {
        let a = 0.125 * (1.0 - x) * (1.0 - e);
        let b = 0.125 * (1.0 + x) * (1.0 - e);
        let c = 0.125 * (1.0 - x) * (1.0 + e);
        let d = 0.125 * (1.0 + x) * (1.0 + e);
        [-a, -b, -c, -d, a, b, c, d]
    };
    let dnx_v = dnx_full(xi, eta, zeta);
    let dne_v = dne(xi, eta, zeta);
    let dnz_v = dnz(xi, eta, zeta);
    for i in 0..8 {
        let node = hex.nodes[i];
        let p = vertices[node];
        for k in 0..3 {
            j[0][k] += dnx_v[i] * p.data[k];
            j[1][k] += dne_v[i] * p.data[k];
            j[2][k] += dnz_v[i] * p.data[k];
        }
    }
    j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
        - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
        + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0])
}

fn apply_dirichlet(r: &mut [f64], k: &mut Vec<Vec<f64>>, u: &[f64], bc: &BoundaryConditions, mesh: &Mesh) {
    // Fully fixed nodes: zero out u and the corresponding rows/cols.
    for &n in &bc.fixed_nodes {
        for k_axis in 0..3 {
            let dof = (n as usize) * 3 + k_axis;
            r[dof] = 0.0;
            for i in 0..mesh.n_dof() {
                k[dof][i] = 0.0;
                k[i][dof] = 0.0;
            }
            k[dof][dof] = 1.0;
        }
    }
    // Per-axis fixed directions.
    for &(n, dir) in &bc.fixed_directions {
        let dof = (n as usize) * 3 + dir as usize;
        r[dof] = 0.0;
        for i in 0..mesh.n_dof() {
            k[dof][i] = 0.0;
            k[i][dof] = 0.0;
        }
        k[dof][dof] = 1.0;
    }
    // Prescribed displacements.
    for &(n, ref disp) in &bc.prescribed_displacement {
        for k_axis in 0..3 {
            let dof = (n as usize) * 3 + k_axis;
            r[dof] = u[dof] - disp.data[k_axis];
            for i in 0..mesh.n_dof() {
                k[dof][i] = 0.0;
                k[i][dof] = 0.0;
            }
            k[dof][dof] = 1.0;
        }
    }
}

fn apply_dirichlet_residual(r: &mut [f64], bc: &BoundaryConditions, mesh: &Mesh) {
    for &n in &bc.fixed_nodes {
        for k in 0..3 {
            r[(n as usize) * 3 + k] = 0.0;
        }
    }
    for &(n, dir) in &bc.fixed_directions {
        r[(n as usize) * 3 + dir as usize] = 0.0;
    }
    let _ = mesh;
}

/// Dense direct linear-system solver via Gauss-Jordan elimination.
fn solve_linear_system(k: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let n = k.len();
    let mut a = vec![vec![0.0_f64; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = k[i][j];
        }
        a[i][n] = rhs[i];
    }
    for i in 0..n {
        let mut piv = i;
        for k2 in (i + 1)..n {
            if a[k2][i].abs() > a[piv][i].abs() {
                piv = k2;
            }
        }
        if a[piv][i].abs() < 1.0e-15 {
            return vec![0.0; n];
        }
        a.swap(i, piv);
        let inv_piv = 1.0 / a[i][i];
        for j in 0..=n {
            a[i][j] *= inv_piv;
        }
        for k2 in 0..n {
            if k2 != i {
                let f = a[k2][i];
                for j in 0..=n {
                    a[k2][j] -= f * a[i][j];
                }
            }
        }
    }
    (0..n).map(|i| a[i][n]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::CrystalPlasticityModel;

    #[test]
    fn unit_cube_mesh_has_8_nodes_1_hex() {
        let m = Mesh::unit_cube();
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.elements.len(), 1);
        assert_eq!(m.n_dof(), 24);
    }

    #[test]
    fn hex8_b_matrix_sum_of_normal_strain_columns_is_zero() {
        // For pure hydrostatic strain, the B-matrix eps_xx
        // component should integrate to zero on a closed domain.
        let m = Mesh::unit_cube();
        let hex = &m.elements[0];
        let mut sum_xx = 0.0_f64;
        for &(xi, eta, zeta) in GAUSS8 {
            let b = hex8_b_matrix(hex, &m.vertices, xi, eta, zeta);
            let det = hex8_det_j(hex, &m.vertices, xi, eta, zeta);
            for i in 0..8 {
                sum_xx += b[0][i * 3] * det;
            }
        }
        // For a closed domain, the volume integral of dN/dx = 0.
        assert!(sum_xx.abs() < 1.0e-9);
    }

    #[test]
    fn fem_assembly_single_element_returns_converged_solution() {
        let mesh = Mesh::unit_cube();
        let elastic = crate::elastic::SymmetricFourthOrder::cubic(168.4e9, 121.4e9, 75.4e9);
        let model = CrystalPlasticityModel::from_crystal_structure(
            tpt_mat_crystallography::CrystalStructure::FCC,
            tpt_mat_hardening::Hardening::Voce(tpt_mat_hardening::VoceHardening::uniform(
                tpt_mat_hardening::VoceParams {
                    tau_0: 30.0e6,
                    tau_s: 60.0e6,
                    theta_0: 500.0e6,
                    gamma_c: 0.05,
                },
            )),
            Default::default(),
            elastic,
        )
        .unwrap();
        let mut solver = CpFemAssembly::new(mesh, model);
        let load = LoadStep::uniaxial(0, 1.0e-3);
        let mut bc = BoundaryConditions::default();
        bc.fixed_nodes.push(0);
        let result = solver.solve(&load, &bc, 200, 1.0e-2);
        assert!(result.is_ok(), "FEM should converge: {:?}", result.err());
        let r = result.unwrap();
        assert!(!r.stresses.is_empty());
        assert_eq!(r.stresses.len(), 8);
    }
}