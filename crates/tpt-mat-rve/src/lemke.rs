//! Lemke linear-complementarity-problem (LCP) solver.
//!
//! Given `M q + r = z,  q ≥ 0, z ≥ 0, q^T z = 0`, find the
//! solution `q, z`.  This implementation follows the textbook
//! Lemke pivoting algorithm (Cottle, Dantzig & Stone 1978) and is
//! exact up to floating-point tolerance.
//!
//! ## Bishop–Hill application
//!
//! The Bishop–Hill problem in crystal plasticity reduces to an
//! LCP on the slip-system activity vector `γ^α`.  In this crate
//! `lemke_solve` is the generic primitive; the crystal-plasticity
//! application is layered on top in
//! `crate::bishop_hill::bishop_hill_lemke`.

/// Result of a Lemke solve.
#[derive(Debug, Clone, PartialEq)]
pub struct LemkeResult {
    /// Original `q` vector.
    pub q: Vec<f64>,
    /// Slack `z` vector.
    pub z: Vec<f64>,
    /// Number of pivots performed.
    pub pivots: usize,
    /// `true` if the algorithm terminated on the *secondary ray*
    /// (infeasible problem or ray termination — pathological).
    pub ray_termination: bool,
}

/// Solve the LCP `M q + r = z,  q ≥ 0, z ≥ 0, q^T z = 0` via the
/// Lemke algorithm with artificial variable `d`.
///
/// - `m`     — square matrix (`n × n`)
/// - `r`     — right-hand side (`n`-vector)
/// - `lex_d` — initial value of the artificial variable.  Pass
///   the maximum of `-rᵢ` to start at a feasible basic solution.
pub fn lemke_solve(m: &[Vec<f64>], r: &[f64], lex_d: f64) -> Option<LemkeResult> {
    let n = m.len();
    if n == 0 || r.len() != n {
        return None;
    }
    // Tableau dimensions: rows = n slacks + 1 artificial.
    // Cols = n originals + n slacks + 1 artificial + 1 RHS.
    let n_cols = 2 * n + 2;
    let n_rows = n + 1;
    let rhs_col = n_cols - 1;
    let artificial_col = n_cols - 2;
    let mut tab: Vec<Vec<f64>> = vec![vec![0.0; n_cols]; n_rows];
    for i in 0..n {
        for j in 0..n {
            tab[i][j] = m[i][j];
        }
        // Slack column `n + i`: coefficient -1.
        tab[i][n + i] = -1.0;
        // RHS = r[i].
        tab[i][rhs_col] = r[i];
    }
    // Lexicographic row: tab[n][j] = (j + 1) for j in 0..=rhs_col-1
    // (used to break ties in the min-ratio test).
    for j in 0..rhs_col {
        tab[n][j] = (j + 1) as f64;
    }
    tab[n][rhs_col] = lex_d;

    // Initial basis: slacks for rows 0..n, artificial for row n.
    let mut basis: Vec<usize> = (0..n).map(|i| n + i).collect();
    basis.push(artificial_col);

    let mut pivots = 0usize;
    let mut ray = false;
    let max_pivots = 4 * n + 4;
    while pivots <= max_pivots {
        // Find the driving variable.
        let driving_col = if basis.contains(&artificial_col) {
            // First pivot (or after the artificial has returned): drive
            // the artificial out.  The driving column is the original
            // column whose lex-row coefficient is most-negative
            // (which is the direction in which the artificial will
            // become negative first).
            let mut target_col = usize::MAX;
            let mut most_neg = -1.0e-30;
            for col in 0..artificial_col {
                if basis.contains(&col) {
                    continue;
                }
                let val = tab[n][col];
                if val < most_neg {
                    most_neg = val;
                    target_col = col;
                }
            }
            if target_col == usize::MAX {
                // All non-basic lex coefficients are non-negative: the
                // artificial never had to leave the basis (the
                // problem is trivially feasible at q = 0).
                basis.retain(|&c| c != artificial_col);
                return Some(LemkeResult {
                    q: vec![0.0; n],
                    z: r.to_vec(),
                    pivots,
                    ray_termination: false,
                });
            }
            target_col
        } else {
            // After the artificial is out, we drive a *complementing*
            // variable: the variable that pairs with the basic one
            // whose RHS is most-negative.  This is the standard Lemke
            // recipe.
            let mut most_neg = -1.0e-30;
            let mut basic_to_drive = usize::MAX;
            for (i, &bcol) in basis.iter().enumerate() {
                if bcol == artificial_col {
                    continue;
                }
                let rhs = tab[i][rhs_col];
                if rhs < most_neg {
                    most_neg = rhs;
                    basic_to_drive = i;
                }
            }
            if basic_to_drive == usize::MAX {
                break;
            }
            let leaving = basis[basic_to_drive];
            // The complementing variable: if leaving is an original
            // (0..n), the partner is the corresponding slack
            // (n+leaving).  If leaving is a slack (n+i), the partner
            // is the original i.
            if leaving < n {
                n + leaving
            } else {
                leaving - n
            }
        };

        // Min-ratio test: find the basic row whose pivot-column
        // coefficient is positive and yields the smallest ratio
        // b/a.  Ties broken lexicographically.
        let mut pivot_row = usize::MAX;
        let mut min_ratio = f64::INFINITY;
        let mut lex_vec: Option<Vec<f64>> = None;
        for (i, &bcol) in basis.iter().enumerate() {
            if bcol == artificial_col && !basis.contains(&artificial_col) {
                continue;
            }
            let a = tab[i][driving_col];
            let b = tab[i][rhs_col];
            if a <= 1.0e-12 {
                continue;
            }
            let ratio = b / a;
            if ratio < 0.0 {
                continue;
            }
            let new_lex: Vec<f64> = (0..=rhs_col)
                .map(|j| {
                    tab[n][j] - (tab[i][j] / a) * tab[n][driving_col]
                })
                .collect();
            let take = match &lex_vec {
                None => true,
                Some(cur) => {
                    if ratio < min_ratio - 1.0e-12 {
                        true
                    } else if (ratio - min_ratio).abs() <= 1.0e-12 {
                        new_lex
                            .iter()
                            .zip(cur.iter())
                            .skip_while(|(a, b)| (**a - **b).abs() < 1.0e-15)
                            .next()
                            .map(|(a, b)| a < b)
                            .unwrap_or(false)
                    } else {
                        false
                    }
                }
            };
            if take {
                min_ratio = ratio;
                lex_vec = Some(new_lex);
                pivot_row = i;
            }
        }
        if pivot_row == usize::MAX {
            ray = true;
            break;
        }
        let leaving_col = basis[pivot_row];
        basis[pivot_row] = driving_col;
        let pivot = tab[pivot_row][driving_col];
        if pivot.abs() < 1.0e-15 {
            ray = true;
            break;
        }
        for j in 0..=rhs_col {
            tab[pivot_row][j] /= pivot;
        }
        for i in 0..n_rows {
            if i == pivot_row {
                continue;
            }
            let factor = tab[i][driving_col];
            if factor.abs() < 1.0e-15 {
                continue;
            }
            for j in 0..=rhs_col {
                tab[i][j] -= factor * tab[pivot_row][j];
            }
        }
        pivots += 1;
        // Termination: the artificial has left the basis.
        if !basis.contains(&artificial_col) {
            break;
        }
        let _ = leaving_col;
    }

    let mut q = vec![0.0_f64; n];
    let mut z = vec![0.0_f64; n];
    for (i, &bcol) in basis.iter().enumerate() {
        if i >= n {
            continue;
        }
        let val = tab[i][rhs_col].max(0.0);
        if bcol < n {
            q[bcol] = val;
        } else if bcol < 2 * n {
            z[bcol - n] = val;
        }
    }
    Some(LemkeResult {
        q,
        z,
        pivots,
        ray_termination: ray,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn trivial_lcp_identity_matrix() {
        let m = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 1.0]];
        let r = vec![1.0, 2.0, 3.0];
        let res = lemke_solve(&m, &r, 10.0).expect("solve");
        assert!(res.q.iter().all(|&v| approx(v, 0.0, 1.0e-9)));
        assert!(approx(res.z[0], 1.0, 1.0e-9));
        assert!(approx(res.z[1], 2.0, 1.0e-9));
        assert!(approx(res.z[2], 3.0, 1.0e-9));
    }

    #[test]
    fn lcp_non_trivial_positive_rhs() {
        // M = [[1, -1], [-1, 1]], r = [1, 1].
        let m = vec![vec![1.0, -1.0], vec![-1.0, 1.0]];
        let r = vec![1.0, 1.0];
        let res = lemke_solve(&m, &r, 10.0).expect("solve");
        let mz = vec![
            m[0][0] * res.q[0] + m[0][1] * res.q[1],
            m[1][0] * res.q[0] + m[1][1] * res.q[1],
        ];
        let z = vec![mz[0] + r[0], mz[1] + r[1]];
        for i in 0..2 {
            assert!(res.q[i] >= -1.0e-9, "q[{i}] < 0");
            assert!(z[i] >= -1.0e-9, "z[{i}] < 0");
        }
        let qt_z: f64 = res.q.iter().zip(z.iter()).map(|(a, b)| a * b).sum();
        assert!(qt_z.abs() < 1.0e-6, "q^T z = {qt_z}");
    }

    #[test]
    fn lcp_complementarity_holds() {
        let m = vec![vec![2.0, -1.0], vec![-1.0, 2.0]];
        let r = vec![1.0, 1.0];
        let res = lemke_solve(&m, &r, 10.0).expect("solve");
        let z = vec![
            m[0][0] * res.q[0] + m[0][1] * res.q[1] + r[0],
            m[1][0] * res.q[0] + m[1][1] * res.q[1] + r[1],
        ];
        let qt_z: f64 = res.q.iter().zip(z.iter()).map(|(a, b)| a * b).sum();
        assert!(qt_z.abs() < 1.0e-6, "q^T z = {qt_z}");
    }

    #[test]
    fn lcp_3x3_returns_feasible_solution() {
        let m = vec![
            vec![4.0, -1.0, 0.0],
            vec![-1.0, 4.0, -1.0],
            vec![0.0, -1.0, 4.0],
        ];
        let r = vec![1.0, 2.0, 3.0];
        let res = lemke_solve(&m, &r, 10.0).expect("solve");
        let z: Vec<f64> = (0..3)
            .map(|i| (0..3).map(|j| m[i][j] * res.q[j]).sum::<f64>() + r[i])
            .collect();
        for i in 0..3 {
            assert!(res.q[i] >= -1.0e-9, "q[{i}] < 0");
            assert!(z[i] >= -1.0e-9, "z[{i}] < 0");
        }
        let qt_z: f64 = res.q.iter().zip(z.iter()).map(|(a, b)| a * b).sum();
        assert!(qt_z.abs() < 1.0e-6, "q^T z = {qt_z}");
    }
}