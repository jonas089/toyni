//! The degree-4 extension `QM31 = CM31[u]/(u^2 - (2 + i))`, the challenge field.
//!
//! A single M31 element gives only ~31 bits of Fiat-Shamir soundness, so all
//! random challenges (constraint batching, the DEEP point, FRI folding
//! randomness) are drawn from this quartic extension of size `p^4 ≈ 2^124`.
//! `2 + i` is a non-residue in CM31, making `u^2 - (2+i)` irreducible; this is
//! the standard "secure field" for M31-based STARKs.
//!
//! Crucially `QM31 ⊇ CM31 = F_p(i)`, so the DEEP single-point quotients of the
//! circle STARK (which are `F(i)`-rational, see Proposition 4 of the paper)
//! live directly in this field — no real/imaginary splitting is needed.

use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use super::cm31::CM31;
use super::m31::M31;

/// `a + b·u` with `u^2 = 2 + i`; coordinates are CM31 elements.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(C)]
pub struct QM31 {
    pub a: CM31,
    pub b: CM31,
}

/// `u^2 = R = 2 + i`.
const R: CM31 = CM31 { a: M31(2), b: M31(1) };

impl QM31 {
    pub const ZERO: Self = Self { a: CM31::ZERO, b: CM31::ZERO };
    pub const ONE: Self = Self { a: CM31::ONE, b: CM31::ZERO };
    /// The imaginary unit `i` of the inner complex extension.
    pub const I: Self = Self { a: CM31::I, b: CM31::ZERO };

    #[inline]
    pub const fn new(a: CM31, b: CM31) -> Self {
        Self { a, b }
    }

    #[inline]
    pub const fn from_cm31(a: CM31) -> Self {
        Self { a, b: CM31::ZERO }
    }

    #[inline]
    pub const fn from_base(a: M31) -> Self {
        Self { a: CM31::from_base(a), b: CM31::ZERO }
    }

    pub fn from_m31_array(c: [M31; 4]) -> Self {
        Self { a: CM31::new(c[0], c[1]), b: CM31::new(c[2], c[3]) }
    }

    pub fn to_m31_array(self) -> [M31; 4] {
        [self.a.a, self.a.b, self.b.a, self.b.b]
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }

    /// True iff the element lies in the base field M31.
    #[inline]
    pub fn is_base(self) -> bool {
        self.a.b.is_zero() && self.b.is_zero()
    }

    #[inline]
    pub fn mul_base(self, s: M31) -> Self {
        Self { a: self.a.mul_base(s), b: self.b.mul_base(s) }
    }

    #[inline]
    pub fn mul_cm31(self, s: CM31) -> Self {
        Self { a: self.a * s, b: self.b * s }
    }

    /// Multiplicative inverse via the quadratic norm down to CM31:
    /// `(a + bu)^-1 = (a - bu) / (a^2 - R·b^2)`.
    pub fn inverse(self) -> Self {
        assert!(!self.is_zero(), "cannot invert zero");
        let denom = self.a * self.a - R * self.b * self.b;
        let denom_inv = denom.inverse();
        Self { a: self.a * denom_inv, b: -self.b * denom_inv }
    }

    pub fn pow(self, mut exp: u128) -> Self {
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

    /// 16-byte little-endian encoding (4 M31 limbs).
    pub fn to_bytes(self) -> [u8; 16] {
        let mut out = [0u8; 16];
        for (i, limb) in self.to_m31_array().iter().enumerate() {
            out[i * 4..(i + 1) * 4].copy_from_slice(&limb.to_bytes());
        }
        out
    }

    pub fn from_bytes(bytes: &[u8; 16]) -> Self {
        let mut c = [M31::ZERO; 4];
        for (i, limb) in c.iter_mut().enumerate() {
            *limb = M31::from_bytes_mod_order(&bytes[i * 4..(i + 1) * 4]);
        }
        Self::from_m31_array(c)
    }

    pub fn random(rng: &mut impl rand::Rng) -> Self {
        Self { a: CM31::random(rng), b: CM31::random(rng) }
    }
}

impl Add for QM31 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        Self { a: self.a + rhs.a, b: self.b + rhs.b }
    }
}

impl Sub for QM31 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        Self { a: self.a - rhs.a, b: self.b - rhs.b }
    }
}

impl Mul for QM31 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        // (a + bu)(c + du) = (ac + R·bd) + (ad + bc)u
        Self {
            a: self.a * rhs.a + R * self.b * rhs.b,
            b: self.a * rhs.b + self.b * rhs.a,
        }
    }
}

impl Neg for QM31 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self { a: -self.a, b: -self.b }
    }
}

impl Div for QM31 {
    type Output = Self;
    #[allow(clippy::suspicious_arithmetic_impl)] // a/b := a·b⁻¹
    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl AddAssign for QM31 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for QM31 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for QM31 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl From<M31> for QM31 {
    fn from(a: M31) -> Self {
        Self::from_base(a)
    }
}

impl From<CM31> for QM31 {
    fn from(a: CM31) -> Self {
        Self::from_cm31(a)
    }
}

impl std::fmt::Display for QM31 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({} + {}u)", self.a, self.b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn rng() -> impl rand::Rng {
        rand::rngs::StdRng::seed_from_u64(11)
    }

    #[test]
    fn u_squared_is_r() {
        let u = QM31::new(CM31::ZERO, CM31::ONE);
        assert_eq!(u * u, QM31::from_cm31(R));
    }

    #[test]
    fn field_axioms_and_inverse() {
        let mut r = rng();
        for _ in 0..500 {
            let a = QM31::random(&mut r);
            let b = QM31::random(&mut r);
            let c = QM31::random(&mut r);
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a * b, b * a);
            assert_eq!((a * b) * c, a * (b * c));
            if !a.is_zero() {
                // Inverses exist for all non-zero elements; this fails if the
                // modulus were reducible (zero divisors are not invertible),
                // so it doubles as an irreducibility check of u^2 - (2+i).
                assert_eq!(a * a.inverse(), QM31::ONE);
            }
        }
    }

    #[test]
    fn tower_embeddings_commute() {
        let mut r = rng();
        for _ in 0..200 {
            let x = M31::random(&mut r);
            let y = M31::random(&mut r);
            assert_eq!(
                QM31::from(x) * QM31::from(y),
                QM31::from(x * y)
            );
            let cx = CM31::random(&mut r);
            let cy = CM31::random(&mut r);
            assert_eq!(QM31::from(cx) * QM31::from(cy), QM31::from(cx * cy));
        }
    }

    #[test]
    fn i_squared_in_qm31() {
        assert_eq!(QM31::I * QM31::I, -QM31::ONE);
    }

    #[test]
    fn byte_roundtrip() {
        let mut r = rng();
        for _ in 0..100 {
            let a = QM31::random(&mut r);
            assert_eq!(QM31::from_bytes(&a.to_bytes()), a);
        }
    }

    #[test]
    fn pow_and_multiplicative_order() {
        // Every nonzero element satisfies x^(p^4 - 1) = 1.
        let p = super::super::m31::P as u128;
        let order = p * p * p * p - 1;
        let mut r = rng();
        for _ in 0..10 {
            let a = QM31::random(&mut r);
            if a.is_zero() {
                continue;
            }
            assert_eq!(a.pow(order), QM31::ONE);
        }
    }
}
