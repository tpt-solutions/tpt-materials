//! Miller indices for crystallographic planes and directions.

use core::fmt;

use serde::{Deserialize, Serialize};

/// Miller indices for a plane or direction.
///
/// A plane is denoted `(h k l)` with the convention that
/// `-1 ≤ h, k, l ≤ 1` for primitive plane families.  A direction is
/// denoted `[u v w]` and may have larger integer components.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct MillerIndex {
    /// First index.
    pub h: i32,
    /// Second index.
    pub k: i32,
    /// Third index.
    pub l: i32,
}

impl MillerIndex {
    /// Construct a Miller index.
    pub const fn new(h: i32, k: i32, l: i32) -> Self {
        Self { h, k, l }
    }
}

impl fmt::Display for MillerIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} {} {})", self.h, self.k, self.l)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        let m = MillerIndex::new(1, 1, 1);
        assert_eq!(format!("{m}"), "(1 1 1)");
    }
}
