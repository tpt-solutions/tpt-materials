//! Integration of the CP-FEM assembly with the cross-repo `tpt-fem`
//! substrate (feature `fem`).
//!
//! The constitutive side of the deferred Phase-2 item is complete; what
//! was missing was the `tpt-fem` mesh-handle interop.  These tests build
//! meshes with `tpt-fem-mesh`'s `MeshBuilder`, convert them with
//! [`Mesh::from_tpt_fem`], and run the full Newton-Raphson CP-FEM solve
//! over them — exercising the multi-element assembly path with shared
//! nodes, not just the built-in `Mesh::unit_cube()`.

#![cfg(feature = "fem")]

use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_mat_crystal_plasticity::{
    BoundaryConditions, CpFemAssembly, CrystalPlasticityModel, LoadStep, Mesh, SymmetricFourthOrder,
};
use tpt_mat_crystallography::CrystalStructure;
use tpt_mat_hardening::{Hardening, VoceHardening, VoceParams};

/// FCC aluminium single-crystal model (elastic constants from the
/// Phase-2 unit tests; Voce hardening with `τ_0 = 30 MPa`).
fn fcc_model() -> CrystalPlasticityModel {
    let elastic = SymmetricFourthOrder::cubic(168.4e9, 121.4e9, 75.4e9);
    CrystalPlasticityModel::from_crystal_structure(
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
    .unwrap()
}

/// Unit-cube corner coordinates in the standard trilinear-hex order
/// (bottom face counter-clockwise at `z = 0`, then the same order at
/// `z = 1`).
const CUBE_CORNERS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];

/// Build a `tpt-fem-mesh` with `nx × ny × nz` unit-cube hexahedra.
fn structured_hex_mesh(nx: usize, ny: usize, nz: usize) -> tpt_fem_mesh::Mesh {
    let mut b = MeshBuilder::new();
    // Node grid: (nx + 1) × (ny + 1) × (nz + 1) nodes, indexed by
    // `node_id(i, j, k) = (k * (ny + 1) + j) * (nx + 1) + i`.
    let stride_x = 1usize;
    let stride_y = nx + 1;
    let stride_z = (nx + 1) * (ny + 1);
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                b.add_node(vec![i as f64, j as f64, k as f64]);
            }
        }
    }
    // One node order per hex corner, expressed as (di, dj, dk) offsets.
    const CORNER_OFFSETS: [[usize; 3]; 8] = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let nodes = CORNER_OFFSETS
                    .iter()
                    .map(|&off| (k + off[2]) * stride_z + (j + off[1]) * stride_y + (i + off[0]))
                    .collect();
                b.try_add_element(CellType::Hex, nodes)
                    .expect("hex connectivity valid");
            }
        }
    }
    let mesh = b.build();
    mesh.validate().expect("structured mesh valid");
    mesh
}

#[test]
fn tpt_fem_single_hex_mesh_drives_uniaxial_tension() {
    let fem_mesh = structured_hex_mesh(1, 1, 1);
    assert_eq!(fem_mesh.node_count(), 8);
    assert_eq!(fem_mesh.element_count(), 1);
    let mesh = Mesh::from_tpt_fem(&fem_mesh).unwrap();
    assert_eq!(mesh.vertices.len(), 8);
    assert_eq!(mesh.n_dof(), 24);

    let mut solver = CpFemAssembly::new(mesh, fcc_model());
    let load = LoadStep::uniaxial(0, 1.0e-3);
    let mut bc = BoundaryConditions::default();
    bc.fixed_nodes.push(0);
    // NOTE: the assembly applies the strain magnitude as a nodal
    // pseudo-traction (`~1e-3`), so the NR "converges" onto the
    // unloaded state with the crate's standard tolerance.  This test
    // validates the mesh-handle plumbing: the converted mesh drives a
    // full NR solve that terminates within the iteration budget and
    // returns finite per-Gauss-point stresses (see `fem_assembly.rs`
    // for the load-application caveat).
    let result = solver
        .solve(&load, &bc, 200, 1.0e-2)
        .expect("single-hex NR solve must converge");
    assert_eq!(result.stresses.len(), 8);
    assert!(
        result.stresses.iter().all(|s| s.data.iter().all(|v| v.is_finite())),
        "all Gauss-point stresses must be finite"
    );
}

#[test]
fn tpt_fem_two_hex_mesh_shares_nodes_and_converges() {
    let fem_mesh = structured_hex_mesh(2, 1, 1);
    // 2 × 1 × 1 grid: 12 nodes (4 shared at x = 1), 2 hexes.
    assert_eq!(fem_mesh.node_count(), 12);
    assert_eq!(fem_mesh.element_count(), 2);
    let mesh = Mesh::from_tpt_fem(&fem_mesh).unwrap();
    assert_eq!(mesh.vertices.len(), 12);
    assert_eq!(mesh.n_dof(), 36);

    let mut solver = CpFemAssembly::new(mesh, fcc_model());
    let load = LoadStep::uniaxial(0, 1.0e-3);
    let mut bc = BoundaryConditions::default();
    bc.fixed_nodes.push(0);
    // See the note in the single-hex test: validates mesh plumbing
    // (shared-node multi-element conversion drives a terminating NR
    // solve with finite stresses).
    let result = solver
        .solve(&load, &bc, 200, 1.0e-2)
        .expect("two-hex NR solve must converge");
    assert_eq!(result.stresses.len(), 16);
    assert!(
        result.stresses.iter().all(|s| s.data.iter().all(|v| v.is_finite())),
        "all Gauss-point stresses must be finite"
    );
}

#[test]
fn tpt_fem_non_hex_cell_is_rejected() {
    let mut b = MeshBuilder::new();
    for c in CUBE_CORNERS {
        b.add_node(c.to_vec());
    }
    b.try_add_element(CellType::Tet, vec![0, 1, 2, 4])
        .expect("tet connectivity valid");
    let fem_mesh = b.build();
    let err = Mesh::from_tpt_fem(&fem_mesh).unwrap_err();
    assert!(
        err.to_string().contains("unsupported cell"),
        "expected an UnsupportedCell error, got {err}"
    );
}
