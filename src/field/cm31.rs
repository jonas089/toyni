//! The complex extension `CM31 = F_p[i]/(i^2 + 1)` of the Mersenne-31 field.
//!
//! Because `p ≡ 3 (mod 4)`, `x^2 + 1` is irreducible over `F_p`, so this is a
//! quadratic field extension of size `p^2`. Geometrically, `i` is the
//! "coordinate at infinity" of the circle curve: the two points at infinity
//! `(1 : ±i : 0)` are CM31-rational, which is what makes single-point (DEEP)
//! quotients work without leaving our challenge field tower.

use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use super::m31::M31;

/// `a + b·i` with `i^2 = -1`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(C)]
pub struct CM31 {
    pub a: M31,
    pub b: M31,
}

impl CM31 {
    pub const ZERO: Self = Self { a: M31::ZERO, b: M31::ZERO };
    pub const ONE: Self = Self { a: M31::ONE, b: M31::ZERO };
    /// The imaginary unit.
    pub const I: Self = Self { a: M31::ZERO, b: M31::ONE };

    #[inline]
    pub const fn new(a: M31, b: M31) -> Self {
        Self { a, b }
    }

    #[inline]
    pub const fn from_base(a: M31) -> Self {
        Self { a, b: M31::ZERO }
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }

    /// Complex conjugate `a - b·i` (the Frobenius `x -> x^p` of CM31/M31).
    #[inline]
    pub fn conjugate(self) -> Self {
        Self { a: self.a, b: -self.b }
    }

    /// The norm `a^2 + b^2 = z · z̄`, an M31 element.
    #[inline]
    pub fn norm(self) -> M31 {
        self.a * self.a + self.b * self.b
    }

    /// Multiplicative inverse: `z^-1 = z̄ / (z·z̄)`.
    pub fn inverse(self) -> Self {
        assert!(!self.is_zero(), "cannot invert zero");
        let norm_inv = self.norm().inverse();
        Self { a: self.a * norm_inv, b: -self.b * norm_inv }
    }

    #[inline]
    pub fn mul_base(self, s: M31) -> Self {
        Self { a: self.a * s, b: self.b * s }
    }

    pub fn random(rng: &mut impl rand::Rng) -> Self {
        Self { a: M31::random(rng), b: M31::random(rng) }
    }
}

impl Add for CM31 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        Self { a: self.a + rhs.a, b: self.b + rhs.b }
    }
}

impl Sub for CM31 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        Self { a: self.a - rhs.a, b: self.b - rhs.b }
    }
}

impl Mul for CM31 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        Self {
            a: self.a * rhs.a - self.b * rhs.b,
            b: self.a * rhs.b + self.b * rhs.a,
        }
    }
}

impl Neg for CM31 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self { a: -self.a, b: -self.b }
    }
}

impl Div for CM31 {
    type Output = Self;
    #[allow(clippy::suspicious_arithmetic_impl)] // a/b := a·b⁻¹
    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl AddAssign for CM31 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for CM31 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for CM31 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl From<M31> for CM31 {
    fn from(a: M31) -> Self {
        Self::from_base(a)
    }
}

impl std::fmt::Display for CM31 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({} + {}i)", self.a, self.b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn rng() -> impl rand::Rng {
        rand::rngs::StdRng::seed_from_u64(7)
    }

    #[test]
    fn i_squared_is_minus_one() {
        assert_eq!(CM31::I * CM31::I, -CM31::ONE);
    }

    #[test]
    fn inverse_and_axioms() {
        let mut r = rng();
        for _ in 0..500 {
            let a = CM31::random(&mut r);
            let b = CM31::random(&mut r);
            let c = CM31::random(&mut r);
            assert_eq!(a * (b + c), a * b + a * c);
            if !a.is_zero() {
                assert_eq!(a * a.inverse(), CM31::ONE);
            }
        }
    }

    #[test]
    fn conjugation_is_ring_hom_and_norm_is_real() {
        let mut r = rng();
        for _ in 0..200 {
            let a = CM31::random(&mut r);
            let b = CM31::random(&mut r);
            assert_eq!((a * b).conjugate(), a.conjugate() * b.conjugate());
            assert_eq!(a * a.conjugate(), CM31::from_base(a.norm()));
        }
    }
}
