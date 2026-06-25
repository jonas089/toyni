//! Degree-4 binomial extension field of BabyBear: `F_p[X] / (X^4 - 11)`.
//!
//! BabyBear is a ~31-bit field, so a single Fiat-Shamir challenge drawn from it
//! only gives ~31 bits of soundness. STARK soundness arguments (lookup/permutation
//! challenges, the out-of-domain point, FRI folding randomness) need much more, so
//! they are drawn from this quartic extension (~124 bits) instead. The trace and
//! the evaluation domain stay in the base field — only random challenges and the
//! values derived from them (accumulators, quotient, DEEP, FRI layers) live here.
//!
//! Element `a0 + a1·X + a2·X^2 + a3·X^3` is stored as `[a0, a1, a2, a3]`. The
//! reduction uses `X^4 = 11` (W = 11; `X^4 - 11` is irreducible over BabyBear,
//! the standard quartic for this field).

use std::hash::{Hash, Hasher};
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::babybear::{BabyBear, BABYBEAR_PRIME};

/// The binomial constant: `X^4 = W`.
const W: u64 = 11;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
pub struct Ext {
    pub c: [BabyBear; 4],
}

impl Hash for Ext {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for limb in &self.c {
            limb.value.hash(state);
        }
    }
}

impl Ext {
    #[inline]
    pub fn new(c: [BabyBear; 4]) -> Self {
        Self { c }
    }

    #[inline]
    pub fn zero() -> Self {
        Self { c: [BabyBear::zero(); 4] }
    }

    #[inline]
    pub fn one() -> Self {
        Self { c: [BabyBear::one(), BabyBear::zero(), BabyBear::zero(), BabyBear::zero()] }
    }

    /// Embed a base-field element as `a + 0·X + 0·X^2 + 0·X^3`.
    #[inline]
    pub fn from_base(b: BabyBear) -> Self {
        Self { c: [b, BabyBear::zero(), BabyBear::zero(), BabyBear::zero()] }
    }

    #[inline]
    pub fn from_u32(value: u32) -> Self {
        Self::from_base(BabyBear::from_u32(value))
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.c.iter().all(|x| x.is_zero())
    }

    /// True iff this element lies in the base field (X-coefficients are zero).
    #[inline]
    pub fn is_base(&self) -> bool {
        self.c[1].is_zero() && self.c[2].is_zero() && self.c[3].is_zero()
    }

    /// Multiply by a base-field scalar (cheaper than a full Ext multiply).
    #[inline]
    pub fn mul_base(self, s: BabyBear) -> Self {
        Self { c: [self.c[0] * s, self.c[1] * s, self.c[2] * s, self.c[3] * s] }
    }

    /// 32-byte little-endian serialization (4 × 8 bytes), for Merkle leaves /
    /// transcript absorption.
    #[inline]
    pub fn to_bytes(&self) -> [u8; 32] {
        let mut out = [0u8; 32];
        for (i, limb) in self.c.iter().enumerate() {
            out[i * 8..i * 8 + 8].copy_from_slice(&limb.to_bytes());
        }
        out
    }

    #[inline]
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut c = [BabyBear::zero(); 4];
        for (i, limb) in c.iter_mut().enumerate() {
            *limb = BabyBear::from_bytes(&bytes[i * 8..i * 8 + 8]);
        }
        Self { c }
    }

    pub fn random(rng: &mut impl rand::Rng) -> Self {
        Self { c: [
            BabyBear::random(rng), BabyBear::random(rng),
            BabyBear::random(rng), BabyBear::random(rng),
        ] }
    }

    /// Square-and-multiply with a u128 exponent (the inverse exponent `p^4 - 2`
    /// exceeds u64).
    pub fn pow_u128(self, mut exp: u128) -> Self {
        let mut base = self;
        let mut result = Self::one();
        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            exp >>= 1;
        }
        result
    }

    /// Multiplicative inverse via Fermat: `a^(p^4 - 2)`.
    pub fn inverse(self) -> Self {
        assert!(!self.is_zero(), "Cannot invert zero");
        let p = BABYBEAR_PRIME as u128;
        let order = p * p * p * p; // p^4
        self.pow_u128(order - 2)
    }
}

impl From<BabyBear> for Ext {
    #[inline]
    fn from(b: BabyBear) -> Self {
        Ext::from_base(b)
    }
}

impl Add for Ext {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self { c: [
            self.c[0] + rhs.c[0], self.c[1] + rhs.c[1],
            self.c[2] + rhs.c[2], self.c[3] + rhs.c[3],
        ] }
    }
}

impl AddAssign for Ext {
    #[inline]
    fn add_assign(&mut self, rhs: Self) { *self = *self + rhs; }
}

impl Sub for Ext {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self { c: [
            self.c[0] - rhs.c[0], self.c[1] - rhs.c[1],
            self.c[2] - rhs.c[2], self.c[3] - rhs.c[3],
        ] }
    }
}

impl SubAssign for Ext {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) { *self = *self - rhs; }
}

impl Neg for Ext {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self { c: [-self.c[0], -self.c[1], -self.c[2], -self.c[3]] }
    }
}

impl Mul for Ext {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let a = &self.c;
        let b = &rhs.c;
        let w = BabyBear::new(W);
        // Schoolbook product mod (X^4 - W), using X^4 = W, X^5 = W·X, X^6 = W·X^2.
        let r0 = a[0] * b[0] + w * (a[1] * b[3] + a[2] * b[2] + a[3] * b[1]);
        let r1 = a[0] * b[1] + a[1] * b[0] + w * (a[2] * b[3] + a[3] * b[2]);
        let r2 = a[0] * b[2] + a[1] * b[1] + a[2] * b[0] + w * (a[3] * b[3]);
        let r3 = a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0];
        Self { c: [r0, r1, r2, r3] }
    }
}

impl MulAssign for Ext {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) { *self = *self * rhs; }
}

impl Div for Ext {
    type Output = Self;
    fn div(self, rhs: Self) -> Self { self * rhs.inverse() }
}

impl std::fmt::Display for Ext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({} + {}x + {}x^2 + {}x^3)", self.c[0], self.c[1], self.c[2], self.c[3])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng() -> impl rand::Rng {
        use rand::SeedableRng;
        rand::rngs::StdRng::seed_from_u64(0xC0FFEE)
    }

    #[test]
    fn base_embedding_is_a_ring_hom() {
        let mut r = rng();
        for _ in 0..200 {
            let x = BabyBear::random(&mut r);
            let y = BabyBear::random(&mut r);
            assert_eq!(Ext::from(x) + Ext::from(y), Ext::from(x + y));
            assert_eq!(Ext::from(x) * Ext::from(y), Ext::from(x * y));
        }
    }

    #[test]
    fn x_to_the_fourth_is_w() {
        // X = [0,1,0,0]; X^4 should equal the base constant W.
        let x = Ext::new([BabyBear::zero(), BabyBear::one(), BabyBear::zero(), BabyBear::zero()]);
        let x4 = x * x * x * x;
        assert_eq!(x4, Ext::from_u32(W as u32));
    }

    #[test]
    fn field_axioms_and_inverse() {
        let mut r = rng();
        // Random elements have inverses (this fails for a reducible modulus,
        // because zero divisors have no inverse — so it also checks irreducibility).
        for _ in 0..500 {
            let a = Ext::random(&mut r);
            if a.is_zero() { continue; }
            assert_eq!(a * a.inverse(), Ext::one());
        }
        // Distributivity.
        for _ in 0..200 {
            let a = Ext::random(&mut r);
            let b = Ext::random(&mut r);
            let c = Ext::random(&mut r);
            assert_eq!(a * (b + c), a * b + a * c);
        }
    }

    #[test]
    fn mul_base_matches_full_mul() {
        let mut r = rng();
        for _ in 0..200 {
            let a = Ext::random(&mut r);
            let s = BabyBear::random(&mut r);
            assert_eq!(a.mul_base(s), a * Ext::from(s));
        }
    }

    #[test]
    fn byte_roundtrip() {
        let mut r = rng();
        for _ in 0..200 {
            let a = Ext::random(&mut r);
            assert_eq!(Ext::from_bytes(&a.to_bytes()), a);
        }
    }
}
