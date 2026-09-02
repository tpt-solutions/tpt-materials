//! Crystal structure enum + lattice parameters.

use serde::{Deserialize, Serialize};

use crate::{LatticeParameters, MillerIndex, SlipSystem};

/// Crystal structure of a phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CrystalStructure {
    /// Face-centred cubic.
    FCC,
    /// Body-centred cubic.
    BCC,
    /// Hexagonal close-packed.
    HCP,
    /// Diamond cubic.
    Diamond,
    /// Simple (primitive) cubic.
    SimpleCubic,
    /// Body-centred tetragonal.
    BCT,
    /// User-supplied structure (no built-in slip systems).
    Custom,
}

impl CrystalStructure {
    /// Default lattice parameters for the structure.  `None` for
    /// `Custom` (the user must provide their own).
    pub fn default_lattice_parameters(self) -> Option<LatticeParameters> {
        // Default values are dimensionless; callers scale with their
        // own lattice parameter (typical FCC Al: 4.05 Å, BCC Fe: 2.87 Å,
        // HCP Ti: a = 2.95 Å, c = 4.69 Å).
        match self {
            CrystalStructure::FCC | CrystalStructure::Diamond | CrystalStructure::SimpleCubic => {
                Some(LatticeParameters::cubic(1.0))
            }
            CrystalStructure::BCC | CrystalStructure::BCT => Some(LatticeParameters::cubic(1.0)),
            CrystalStructure::HCP => Some(LatticeParameters::hcp(1.0, 1.633)),
            CrystalStructure::Custom => None,
        }
    }

    /// Slip systems for the structure.
    ///
    /// | Structure | Count | Family |
    /// |---|---|---|
    /// | FCC | 12 | `{111}⟨110⟩` |
    /// | BCC | 12 | `{110}⟨111⟩` |
    /// | HCP | 3 basal + 3 prismatic + 6 pyramidal ⟨a⟩ + 12 pyramidal ⟨c+a⟩ = 24 |
    /// | Diamond | 0 (use FCC slip systems by convention) |
    /// | SimpleCubic | 0 |
    /// | BCT | 0 |
    /// | Custom | 0 |
    pub fn slip_systems(self) -> Vec<SlipSystem> {
        match self {
            CrystalStructure::FCC => fcc_slip_systems(),
            CrystalStructure::BCC => bcc_slip_systems(),
            CrystalStructure::HCP => hcp_slip_systems(),
            _ => Vec::new(),
        }
    }
}

/// FCC: 12 slip systems, family `{111}⟨110⟩`.
fn fcc_slip_systems() -> Vec<SlipSystem> {
    // 4 planes × 3 directions per plane = 12.
    let planes: [MillerIndex; 4] = [
        MillerIndex::new(1, 1, 1),
        MillerIndex::new(-1, 1, 1),
        MillerIndex::new(1, -1, 1),
        MillerIndex::new(1, 1, -1),
    ];
    let dirs_per_plane: [[i32; 3]; 3] = [[-1, 1, 0], [0, -1, 1], [1, 0, -1]];
    // Combine, then dedupe by (plane, direction) up to sign symmetry.
    let mut systems = Vec::with_capacity(12);
    for plane in &planes {
        for d in &dirs_per_plane {
            let n = plane_normal(*plane);
            let s = slip_direction([d[0], d[1], d[2]]);
            systems.push(SlipSystem {
                plane: *plane,
                slip_direction_miller: *d,
                plane_normal: n,
                slip_direction: s,
                critical_resolved_shear_stress: 1.0,
                family: crate::SlipFamily::Fcc110,
            });
        }
    }
    systems
}

/// BCC: 12 slip systems, family `{110}⟨111⟩`.
fn bcc_slip_systems() -> Vec<SlipSystem> {
    let planes: [MillerIndex; 6] = [
        MillerIndex::new(1, 1, 0),
        MillerIndex::new(-1, 1, 0),
        MillerIndex::new(1, 0, 1),
        MillerIndex::new(1, 0, -1),
        MillerIndex::new(0, 1, 1),
        MillerIndex::new(0, -1, 1),
    ];
    // 2 unique ⟨111⟩ directions per plane (up to sign).
    let dirs_per_plane: [[i32; 3]; 2] = [[1, 1, 1], [-1, -1, 1]];
    let mut systems = Vec::with_capacity(12);
    for plane in &planes {
        for d in &dirs_per_plane {
            let n = plane_normal(*plane);
            let s = slip_direction([d[0], d[1], d[2]]);
            systems.push(SlipSystem {
                plane: *plane,
                slip_direction_miller: *d,
                plane_normal: n,
                slip_direction: s,
                critical_resolved_shear_stress: 1.0,
                family: crate::SlipFamily::Bcc111,
            });
        }
    }
    systems
}

/// HCP: basal ⟨a⟩ + prismatic ⟨a⟩ + pyramidal ⟨a⟩ + pyramidal ⟨c+a⟩.
fn hcp_slip_systems() -> Vec<SlipSystem> {
    let mut systems = Vec::with_capacity(24);

    // Basal: 3 systems, family (0001)<11-20>.
    let basal_planes = [
        MillerIndex::new(0, 0, 0),
        MillerIndex::new(0, 0, 0),
        MillerIndex::new(0, 0, 0),
    ];
    let basal_dirs: [[i32; 3]; 3] = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]];
    for (i, d) in basal_dirs.iter().enumerate() {
        systems.push(SlipSystem {
            plane: basal_planes[i],
            slip_direction_miller: *d,
            // For HCP basal we represent the plane normal by the c-axis
            // (Miller-Bravais [0001] → 3-index [001]) and the
            // <11-20> directions as the closest 3-index <a> vectors.
            plane_normal: plane_normal_3index([0, 0, 1]),
            slip_direction: slip_direction([d[0], d[1], d[2]]),
            critical_resolved_shear_stress: 1.0,
            family: crate::SlipFamily::HcpBasal,
        });
    }

    // Prismatic: 3 systems, family {10-10}<11-20>.
    let pris_planes: [MillerIndex; 3] = [
        MillerIndex::new(1, -1, 0),
        MillerIndex::new(0, 1, -1),
        MillerIndex::new(-1, 0, 1),
    ];
    for (i, d) in basal_dirs.iter().enumerate() {
        let p = pris_planes[i];
        systems.push(SlipSystem {
            plane: p,
            slip_direction_miller: *d,
            plane_normal: plane_normal_3index([p.h, p.k, p.l]),
            slip_direction: slip_direction([d[0], d[1], d[2]]),
            critical_resolved_shear_stress: 1.0,
            family: crate::SlipFamily::HcpPrismatic,
        });
    }

    // Pyramidal ⟨a⟩: 6 systems, family {10-11}<11-20>.
    let pyr_a_planes: [MillerIndex; 6] = [
        MillerIndex::new(1, -1, 1),
        MillerIndex::new(-1, 1, 1),
        MillerIndex::new(1, 0, 1),
        MillerIndex::new(-1, 0, 1),
        MillerIndex::new(0, 1, 1),
        MillerIndex::new(0, -1, 1),
    ];
    for (i, p) in pyr_a_planes.iter().enumerate() {
        let d = basal_dirs[i % 3];
        systems.push(SlipSystem {
            plane: *p,
            slip_direction_miller: d,
            plane_normal: plane_normal_3index([p.h, p.k, p.l]),
            slip_direction: slip_direction([d[0], d[1], d[2]]),
            critical_resolved_shear_stress: 1.0,
            family: crate::SlipFamily::HcpPyramidalA,
        });
    }

    // Pyramidal ⟨c+a⟩: 12 systems, family {11-22}<11-23>.
    let pyr_ca_planes: [MillerIndex; 12] = [
        MillerIndex::new(1, 1, -2),
        MillerIndex::new(-1, -1, -2),
        MillerIndex::new(1, -1, -2),
        MillerIndex::new(-1, 1, -2),
        MillerIndex::new(1, 0, -2),
        MillerIndex::new(-1, 0, -2),
        MillerIndex::new(0, 1, -2),
        MillerIndex::new(0, -1, -2),
        MillerIndex::new(1, 1, 2),
        MillerIndex::new(-1, -1, 2),
        MillerIndex::new(1, -1, 2),
        MillerIndex::new(-1, 1, 2),
    ];
    let ca_dirs: [[i32; 3]; 12] = [
        [-2, 1, 1],
        [2, -1, -1],
        [1, -2, 1],
        [-1, 2, -1],
        [-2, 0, 1],
        [2, 0, -1],
        [0, -2, 1],
        [0, 2, -1],
        [1, 1, -2],
        [-1, -1, 2],
        [2, -1, -2],
        [-2, 1, 2],
    ];
    for (i, p) in pyr_ca_planes.iter().enumerate() {
        let d = ca_dirs[i];
        systems.push(SlipSystem {
            plane: *p,
            slip_direction_miller: d,
            plane_normal: plane_normal_3index([p.h, p.k, p.l]),
            slip_direction: slip_direction([d[0], d[1], d[2]]),
            critical_resolved_shear_stress: 1.0,
            family: crate::SlipFamily::HcpPyramidalCA,
        });
    }

    systems
}

fn plane_normal(m: MillerIndex) -> tpt_math_linalg_fixed::Vec3 {
    plane_normal_3index([m.h, m.k, m.l])
}

fn plane_normal_3index(m: [i32; 3]) -> tpt_math_linalg_fixed::Vec3 {
    tpt_math_linalg_fixed::Vec3::new(m[0] as f64, m[1] as f64, m[2] as f64).normalized()
}

fn slip_direction(d: [i32; 3]) -> tpt_math_linalg_fixed::Vec3 {
    tpt_math_linalg_fixed::Vec3::new(d[0] as f64, d[1] as f64, d[2] as f64).normalized()
}
