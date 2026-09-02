//! Regular Cartesian grids + finite-difference / spectral Laplacian.
//!
//! Used by `tpt-mat-phase-field` (Allen-Cahn / Cahn-Hilliard),
//! `tpt-mat-grain-growth`, `tpt-mat-solidification`, and
//! `tpt-mat-diffusion`.  All grids are row-major with `dx` the uniform
//! cell spacing.

use serde::{Deserialize, Serialize};

/// 1D regular grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Grid1D {
    /// Number of cells.
    pub nx: usize,
    /// Cell spacing.
    pub dx: f64,
}

impl Grid1D {
    /// Construct a 1D grid.
    pub fn new(nx: usize, dx: f64) -> Self {
        Self { nx, dx }
    }

    /// Length of the domain.
    pub fn length(&self) -> f64 {
        self.nx as f64 * self.dx
    }

    /// Apply the standard 3-point Laplacian to a 1D field `u`.  The
    /// boundary values are passed in (Dirichlet) and left untouched.
    pub fn laplacian(&self, u: &[f64], out: &mut [f64], left: f64, right: f64) {
        assert_eq!(u.len(), self.nx);
        assert_eq!(out.len(), self.nx);
        let inv_dx2 = 1.0 / (self.dx * self.dx);
        for i in 0..self.nx {
            let ul = if i == 0 { left } else { u[i - 1] };
            let ur = if i + 1 == self.nx { right } else { u[i + 1] };
            out[i] = (ul - 2.0 * u[i] + ur) * inv_dx2;
        }
    }

    /// Apply the biharmonic operator `∇⁴ = ∇²∇²` (used by the
    /// Cahn-Hilliard equation) to a 1D field.  Boundaries are
    /// zero-flux (Neumann).
    pub fn biharmonic(&self, u: &[f64], out: &mut [f64]) {
        let mut lap = vec![0.0_f64; self.nx];
        // Neumann: use one-sided stencil at the boundaries.
        self.laplacian_neumann(u, &mut lap);
        self.laplacian_neumann(&lap, out);
    }

    /// Laplacian with zero-flux Neumann boundaries.
    pub fn laplacian_neumann(&self, u: &[f64], out: &mut [f64]) {
        assert_eq!(u.len(), self.nx);
        assert_eq!(out.len(), self.nx);
        let inv_dx2 = 1.0 / (self.dx * self.dx);
        for i in 0..self.nx {
            let ul = if i == 0 { u[i + 1] } else { u[i - 1] };
            let ur = if i + 1 == self.nx { u[i - 1] } else { u[i + 1] };
            out[i] = (ul - 2.0 * u[i] + ur) * inv_dx2;
        }
    }
}

/// 2D regular grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Grid2D {
    /// Cells along x.
    pub nx: usize,
    /// Cells along y.
    pub ny: usize,
    /// Cell spacing (square cells assumed).
    pub dx: f64,
}

impl Grid2D {
    /// Construct a square 2D grid.
    pub fn new(nx: usize, ny: usize, dx: f64) -> Self {
        Self { nx, ny, dx }
    }

    /// Total number of cells.
    pub fn len(&self) -> usize {
        self.nx * self.ny
    }

    /// Whether the grid is empty.
    pub fn is_empty(&self) -> bool {
        self.nx == 0 || self.ny == 0
    }

    /// Convert `(i, j)` to a flat index (row-major, `i` is row, `j` is column).
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        i * self.nx + j
    }

    /// Apply the 5-point Laplacian with zero-flux Neumann boundaries.
    /// Inputs / outputs are flat length-`nx*ny` row-major arrays.
    pub fn laplacian_neumann(&self, u: &[f64], out: &mut [f64]) {
        let n = self.len();
        assert_eq!(u.len(), n);
        assert_eq!(out.len(), n);
        let inv_dx2 = 1.0 / (self.dx * self.dx);
        for i in 0..self.ny {
            for j in 0..self.nx {
                let c = self.idx(i, j);
                let u_c = u[c];
                let up = if i + 1 < self.ny { u[self.idx(i + 1, j)] } else { u[self.idx(i - 1, j)] };
                let um = if i > 0 { u[self.idx(i - 1, j)] } else { u[self.idx(i + 1, j)] };
                let ur = if j + 1 < self.nx { u[self.idx(i, j + 1)] } else { u[self.idx(i, j - 1)] };
                let ul = if j > 0 { u[self.idx(i, j - 1)] } else { u[self.idx(i, j + 1)] };
                out[c] = (up + um + ur + ul - 4.0 * u_c) * inv_dx2;
            }
        }
    }

    /// Apply the biharmonic operator with zero-flux Neumann boundaries.
    pub fn biharmonic(&self, u: &[f64], out: &mut [f64]) {
        let mut lap = vec![0.0_f64; self.len()];
        self.laplacian_neumann(u, &mut lap);
        self.laplacian_neumann(&lap, out);
    }

    /// Apply the biharmonic operator with periodic boundaries.
    pub fn biharmonic_periodic(&self, u: &[f64], out: &mut [f64]) {
        let mut lap = vec![0.0_f64; self.len()];
        self.laplacian_periodic(u, &mut lap);
        self.laplacian_periodic(&lap, out);
    }

    /// Apply the Laplacian with periodic boundaries.
    pub fn laplacian_periodic(&self, u: &[f64], out: &mut [f64]) {
        let n = self.len();
        assert_eq!(u.len(), n);
        assert_eq!(out.len(), n);
        let inv_dx2 = 1.0 / (self.dx * self.dx);
        for i in 0..self.ny {
            let ip = (i + 1) % self.ny;
            let im = (i + self.ny - 1) % self.ny;
            for j in 0..self.nx {
                let jp = (j + 1) % self.nx;
                let jm = (j + self.nx - 1) % self.nx;
                let c = self.idx(i, j);
                out[c] = (u[self.idx(ip, j)] + u[self.idx(im, j)]
                    + u[self.idx(i, jp)] + u[self.idx(i, jm)]
                    - 4.0 * u[c])
                    * inv_dx2;
            }
        }
    }
}

/// 3D regular grid (used by `tpt-mat-solidification`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Grid3D {
    /// Cells along x.
    pub nx: usize,
    /// Cells along y.
    pub ny: usize,
    /// Cells along z.
    pub nz: usize,
    /// Cell spacing.
    pub dx: f64,
}

impl Grid3D {
    /// Construct a cubic 3D grid.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        Self { nx, ny, nz, dx }
    }

    /// Total number of cells.
    pub fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Whether the grid is empty.
    pub fn is_empty(&self) -> bool {
        self.nx == 0 || self.ny == 0 || self.nz == 0
    }

    /// Flat index `(i, j, k)` in `(z, y, x)` row-major order.
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.ny + j) * self.nx + k
    }

    /// Apply the 7-point Laplacian with zero-flux Neumann boundaries.
    pub fn laplacian_neumann(&self, u: &[f64], out: &mut [f64]) {
        let n = self.len();
        assert_eq!(u.len(), n);
        assert_eq!(out.len(), n);
        let inv_dx2 = 1.0 / (self.dx * self.dx);
        for i in 0..self.nz {
            for j in 0..self.ny {
                for k in 0..self.nx {
                    let c = self.idx(i, j, k);
                    let u_c = u[c];
                    let u_ip = if i + 1 < self.nz {
                        u[self.idx(i + 1, j, k)]
                    } else {
                        u[self.idx(i - 1, j, k)]
                    };
                    let u_im = if i > 0 {
                        u[self.idx(i - 1, j, k)]
                    } else {
                        u[self.idx(i + 1, j, k)]
                    };
                    let u_jp = if j + 1 < self.ny {
                        u[self.idx(i, j + 1, k)]
                    } else {
                        u[self.idx(i, j - 1, k)]
                    };
                    let u_jm = if j > 0 {
                        u[self.idx(i, j - 1, k)]
                    } else {
                        u[self.idx(i, j + 1, k)]
                    };
                    let u_kp = if k + 1 < self.nx {
                        u[self.idx(i, j, k + 1)]
                    } else {
                        u[self.idx(i, j, k - 1)]
                    };
                    let u_km = if k > 0 {
                        u[self.idx(i, j, k - 1)]
                    } else {
                        u[self.idx(i, j, k + 1)]
                    };
                    out[c] = (u_ip + u_im + u_jp + u_jm + u_kp + u_km - 6.0 * u_c) * inv_dx2;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn grid1d_laplacian_of_quadratic() {
        // u(x) = x^2 → u''(x) = 2 on interior cells.
        let g = Grid1D::new(5, 0.1);
        let u: Vec<f64> = (0..5).map(|i| (i as f64 * 0.1).powi(2)).collect();
        let mut lap = vec![0.0_f64; 5];
        g.laplacian_neumann(&u, &mut lap);
        for i in 1..4 {
            assert_relative_eq!(lap[i], 2.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn grid2d_laplacian_periodic_of_sine() {
        // u(x,y) = sin(2π x/L + 2π y/L) → ∇² u = -k² u, k² = 2 (2π/L)².
        let g = Grid2D::new(16, 16, 1.0 / 16.0);
        let l = g.nx as f64 * g.dx;
        let k = 2.0 * std::f64::consts::PI / l;
        let u: Vec<f64> = (0..g.len())
            .map(|c| {
                let i = c / g.nx;
                let j = c % g.nx;
                (k * (i as f64 * g.dx + j as f64 * g.dx)).sin()
            })
            .collect();
        let mut lap = vec![0.0_f64; g.len()];
        g.laplacian_periodic(&u, &mut lap);
        let ksq = 2.0 * k * k;
        // 2nd-order central-difference error scales as O((k·dx)²); with
        // 16 cells per wavelength the relative error is ~1.3%.
        // Compare per-cell absolute error against the analytical max.
        let analytical_max = ksq;
        for (u_c, lap_c) in u.iter().zip(lap.iter()) {
            let expected = -ksq * u_c;
            assert!(
                (*lap_c - expected).abs() / analytical_max < 0.02,
                "lap={lap_c}, expected={expected}, u={u_c}"
            );
        }
    }

    #[test]
    fn grid3d_laplacian_neumann_constant_field_is_zero() {
        let g = Grid3D::new(4, 4, 4, 0.5);
        let u = vec![1.0_f64; g.len()];
        let mut lap = vec![0.0_f64; g.len()];
        g.laplacian_neumann(&u, &mut lap);
        for v in lap {
            assert_eq!(v, 0.0);
        }
    }
}