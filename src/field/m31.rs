//! The Mersenne-31 prime field, `F_p` with `p = 2^31 - 1`.
//!
//! M31 is the canonical CFFT-friendly prime: `p + 1 = 2^31`, so the circle
//! group `C(F_p)` is cyclic of order `2^31` and supports circle-FFT domains
//! of every two-adic size up to `2^30`. Arithmetic is exceptionally fast
//! because reduction modulo `2^31 - 1` is a shift and an add.

use std::hash::{Hash, Hasher};
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// The Mersenne prime `2^31 - 1`.
pub const P: u32 = (1 << 31) - 1;

/// An element of the Mersenne-31 prime field, kept reduced in `[0, P)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(transparent)]
pub struct M31(pub u32);

impl Hash for M31 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl M31 {
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);
    pub const TWO: Self = Self(2);

    /// Reduce an arbitrary `u32` into the field. Note `P` itself maps to 0.
    #[inline(always)]
    pub const fn new(value: u32) -> Self {
        Self(value % P)
    }

    /// Reduce a `u64` (e.g. a product of two reduced elements) into the field.
    ///
    /// Splits into 31-bit limbs: `x = hi * 2^31 + lo ≡ hi + lo (mod p)`.
    #[inline(always)]
    pub const fn reduce_u64(x: u64) -> Self {
        let first = (x >> 31) + (x & P as u64); // < 2^33
        let second = (first >> 31) + (first & P as u64); // <= 2^31
        let mut r = second as u32;
        if r >= P {
            r -= P;
        }
        Self(r)
    }

    #[inline(always)]
    pub const fn is_zero(self) -> bool {
        self.0 == 0
    }

    #[inline]
    pub fn from_u64(value: u64) -> Self {
        Self::reduce_u64(value)
    }

    /// Interpret up to 4 little-endian bytes, reduced mod `p`.
    pub fn from_bytes_mod_order(bytes: &[u8]) -> Self {
        let mut arr = [0u8; 4];
        let n = bytes.len().min(4);
        arr[..n].copy_from_slice(&bytes[..n]);
        Self::new(u32::from_le_bytes(arr))
    }

    #[inline]
    pub fn to_bytes(self) -> [u8; 4] {
        self.0.to_le_bytes()
    }

    pub fn random(rng: &mut impl rand::Rng) -> Self {
        Self(rng.gen_range(0..P))
    }

    /// Exponentiation by squaring.
    pub fn pow(self, mut exp: u64) -> Self {
        let mut base = self;
        let mut result = Self::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                result *= base;
            }
            base *= base;
            exp >>= 1;
        }
        result
    }

    /// Multiplicative inverse via Fermat: `a^(p-2)`.
    ///
    /// Panics on zero.
    pub fn inverse(self) -> Self {
        assert!(!self.is_zero(), "cannot invert zero");
        self.pow(P as u64 - 2)
    }

    /// `1/2` in the field, useful for FFT-style averaging.
    #[inline]
    pub const fn half() -> Self {
        // 2 * 2^30 = 2^31 ≡ 1 (mod p)
        Self(1 << 30)
    }
}

impl Add for M31 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        let mut s = self.0 + rhs.0; // < 2^32, no overflow
        if s >= P {
            s -= P;
        }
        Self(s)
    }
}

impl AddAssign for M31 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for M31 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        let (d, borrow) = self.0.overflowing_sub(rhs.0);
        Self(if borrow { d.wrapping_add(P) } else { d })
    }
}

impl SubAssign for M31 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for M31 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        Self::reduce_u64(self.0 as u64 * rhs.0 as u64)
    }
}

impl MulAssign for M31 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Div for M31 {
    type Output = Self;
    #[allow(clippy::suspicious_arithmetic_impl)] // a/b := a·b⁻¹
    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl Neg for M31 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        if self.0 == 0 {
            self
        } else {
            Self(P - self.0)
        }
    }
}

impl From<u32> for M31 {
    fn from(value: u32) -> Self {
        Self::new(value)
    }
}

impl std::fmt::Display for M31 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn rng() -> impl rand::Rng {
        rand::rngs::StdRng::seed_from_u64(0xC13C1E)
    }

    #[test]
    fn reduction_edge_cases() {
        assert_eq!(M31::new(P), M31::ZERO);
        assert_eq!(M31::new(P - 1) + M31::ONE, M31::ZERO);
        assert_eq!(M31::reduce_u64(u64::MAX), M31::new((u64::MAX % P as u64) as u32));
        assert_eq!(M31::reduce_u64(P as u64 * P as u64), M31::ZERO);
        // (p-1)^2 mod p = 1
        assert_eq!(M31::new(P - 1) * M31::new(P - 1), M31::ONE);
    }

    #[test]
    fn field_axioms_random() {
        let mut r = rng();
        for _ in 0..1000 {
            let a = M31::random(&mut r);
            let b = M31::random(&mut r);
            let c = M31::random(&mut r);
            assert_eq!(a + b, b + a);
            assert_eq!(a * b, b * a);
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a - b + b, a);
            if !a.is_zero() {
                assert_eq!(a * a.inverse(), M31::ONE);
            }
        }
    }

    #[test]
    fn half_is_inverse_of_two() {
        assert_eq!(M31::half() * M31::TWO, M31::ONE);
    }

    #[test]
    fn pow_matches_naive() {
        let a = M31::new(3);
        let mut acc = M31::ONE;
        for e in 0..40u64 {
            assert_eq!(a.pow(e), acc);
            acc *= a;
        }
    }

    #[test]
    fn neg_and_sub() {
        let mut r = rng();
        for _ in 0..100 {
            let a = M31::random(&mut r);
            assert_eq!(a + (-a), M31::ZERO);
            assert_eq!(M31::ZERO - a, -a);
        }
    }
}
