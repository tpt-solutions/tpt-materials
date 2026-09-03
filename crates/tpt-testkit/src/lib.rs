//! Minimal floating-point comparison assertions for workspace tests.
//!
//! This is a small, from-scratch, MIT-licensed replacement for the
//! `approx` crate (which is Apache-2.0-only and would break the pure-MIT
//! dependency path this workspace guarantees — see `LICENSING.md`).
//!
//! It implements only the subset used by the workspace:
//!
//! - [`relative_eq`] / [`abs_diff_eq`] — boolean comparisons
//! - [`assert_relative_eq!`] / [`assert_abs_diff_eq!`] — assertions
//!
//! # Semantics
//!
//! [`relative_eq`] mirrors `approx::relative_eq`: two values compare
//! equal when either
//!
//! ```text
//! |a − b| <= epsilon                              (absolute test)
//! |a − b| <= max_relative * max(|a|, |b|)         (relative test)
//! ```
//!
//! holds.  Both tolerances default to [`f64::EPSILON`].  The macros
//! accept the same `epsilon = …` and `max_relative = …` keyword
//! arguments as `approx`.

/// Boolean absolute-difference comparison: `|a − b| <= epsilon`.
#[must_use]
pub fn abs_diff_eq(a: f64, b: f64, epsilon: f64) -> bool {
    if a == b {
        return true;
    }
    if a.is_infinite() || b.is_infinite() || a.is_nan() || b.is_nan() {
        return false;
    }
    (a - b).abs() <= epsilon
}

/// Boolean relative comparison mirroring `approx::relative_eq`.
///
/// Returns `true` when the absolute test (`|a − b| <= epsilon`) *or* the
/// relative test (`|a − b| <= max_relative * max(|a|, |b|)`) passes.
#[must_use]
pub fn relative_eq(a: f64, b: f64, epsilon: f64, max_relative: f64) -> bool {
    if a == b {
        return true;
    }
    if a.is_infinite() || b.is_infinite() || a.is_nan() || b.is_nan() {
        return false;
    }
    let diff = (a - b).abs();
    if diff <= epsilon {
        return true;
    }
    let largest = a.abs().max(b.abs());
    diff <= largest * max_relative
}

/// Default tolerance used when no keyword argument is supplied.
pub const DEFAULT_EPSILON: f64 = f64::EPSILON;

/// Assert that two `f64` expressions are relatively equal.
///
/// ```
/// # use tpt_testkit::assert_relative_eq;
/// assert_relative_eq!(1.0_f64, 1.0 + 1e-12, epsilon = 1e-9);
/// assert_relative_eq!(1.0e6_f64, 1.0e6 + 1.0, max_relative = 1e-3);
/// ```
#[macro_export]
macro_rules! assert_relative_eq {
    ($a:expr, $b:expr $(,)?) => {
        $crate::assert_relative_eq!(
            $a,
            $b,
            epsilon = $crate::DEFAULT_EPSILON,
            max_relative = $crate::DEFAULT_EPSILON
        )
    };
    ($a:expr, $b:expr, epsilon = $eps:expr $(,)?) => {
        $crate::assert_relative_eq!(
            $a,
            $b,
            epsilon = $eps,
            max_relative = $crate::DEFAULT_EPSILON
        )
    };
    ($a:expr, $b:expr, max_relative = $mr:expr $(,)?) => {
        $crate::assert_relative_eq!(
            $a,
            $b,
            epsilon = $crate::DEFAULT_EPSILON,
            max_relative = $mr
        )
    };
    ($a:expr, $b:expr, epsilon = $eps:expr, max_relative = $mr:expr $(,)?) => {{
        let a: f64 = $a as f64;
        let b: f64 = $b as f64;
        assert!(
            $crate::relative_eq(a, b, $eps, $mr),
            "assert_relative_eq!({}, {})\n     left  = {:?}\n     right = {:?}\n     |diff| = {:?} (epsilon = {:?}, max_relative = {:?})",
            stringify!($a),
            stringify!($b),
            a,
            b,
            (a - b).abs(),
            $eps as f64,
            $mr as f64,
        );
    }};
    ($a:expr, $b:expr, max_relative = $mr:expr, epsilon = $eps:expr $(,)?) => {
        $crate::assert_relative_eq!($a, $b, epsilon = $eps, max_relative = $mr)
    };
}

/// Assert that two `f64` expressions differ by at most `epsilon`
/// (absolute only). Defaults to [`DEFAULT_EPSILON`].
///
/// ```
/// # use tpt_testkit::assert_abs_diff_eq;
/// assert_abs_diff_eq!(0.1_f64 + 0.2, 0.3, epsilon = 1e-12);
/// ```
#[macro_export]
macro_rules! assert_abs_diff_eq {
    ($a:expr, $b:expr $(,)?) => {
        $crate::assert_abs_diff_eq!($a, $b, epsilon = $crate::DEFAULT_EPSILON)
    };
    ($a:expr, $b:expr, epsilon = $eps:expr $(,)?) => {{
        let a: f64 = $a as f64;
        let b: f64 = $b as f64;
        assert!(
            $crate::abs_diff_eq(a, b, $eps),
            "assert_abs_diff_eq!({}, {})\n     left  = {:?}\n     right = {:?}\n     |diff| = {:?} (epsilon = {:?})",
            stringify!($a),
            stringify!($b),
            a,
            b,
            (a - b).abs(),
            $eps as f64,
        );
    }};
}

#[cfg(test)]
mod tests {
    #[test]
    fn abs_and_relative_defaults() {
        assert!(super::abs_diff_eq(1.0, 1.0, 0.0));
        assert!(!super::abs_diff_eq(1.0, 1.1, 1e-3));
        assert!(super::relative_eq(1.0e9, 1.0e9 + 1.0, f64::EPSILON, 1e-6));
        assert!(!super::relative_eq(1.0, 2.0, f64::EPSILON, f64::EPSILON));
    }

    #[test]
    fn macro_forms() {
        assert_relative_eq!(1.0_f64, 1.0);
        assert_relative_eq!(1.0_f64, 1.0 + 1e-12, epsilon = 1e-9);
        assert_relative_eq!(1.0e6_f64, 1.0e6 + 1.0, max_relative = 1e-3);
        assert_relative_eq!(2.0_f64, 2.0 + 1e-12, epsilon = 1e-9, max_relative = 1e-9);
        assert_abs_diff_eq!(0.1_f64 + 0.2, 0.3, epsilon = 1e-12);
    }

    #[test]
    #[should_panic(expected = "assert_relative_eq!")]
    fn macro_fails_loudly() {
        assert_relative_eq!(1.0_f64, 2.0, epsilon = 1e-9);
    }
}
