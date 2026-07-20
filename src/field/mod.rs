//! The field abstraction shared by both proving engines.
//!
//! toyni now hosts two STARK engines over two different fields:
//!
//! - the **classical** roots-of-unity engine over BabyBear
//!   (`p = 2^31 − 2^27 + 1`), with its quartic extension [`babybear_ext::Ext`];
//! - the **circle** engine over Mersenne-31 (`p = 2^31 − 1`), with the field
//!   tower `M31 ⊂ CM31 ⊂ QM31` (the `circle` feature).
//!
//! Both are exposed through one interface so that AIRs, examples, Merkle
//! commitments and the Fiat–Shamir transcript are written once and run on
//! either engine:
//!
//! - [`Field`] — a prime or extension field element (the operations every
//!   committed value and challenge needs);
//! - [`ExtField`] — an extension `E` over base `B` with a cheap base-scalar
//!   multiply and polynomial evaluation, home of the ~124-bit challenges;
//! - [`ConstraintField`] — a field in which an AIR whose constants live in
//!   base `B` can be evaluated: either `B` itself (the LDE hot path) or an
//!   extension of it (the out-of-domain check).

use rand::Rng;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub mod babybear;
pub mod babybear_ext;

#[cfg(feature = "circle")]
pub mod cm31;
#[cfg(feature = "circle")]
pub mod m31;
#[cfg(feature = "circle")]
pub mod qm31;

pub use babybear::BabyBear;
pub use babybear_ext::Ext as BabyBearExt;

#[cfg(feature = "circle")]
pub use cm31::CM31;
#[cfg(feature = "circle")]
pub use m31::M31;
#[cfg(feature = "circle")]
pub use qm31::QM31;

/// A finite-field element usable as a base or extension field in the STARK.
pub trait Field:
    Copy
    + Clone
    + PartialEq
    + Eq
    + Debug
    + Default
    + Hash
    + Send
    + Sync
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
{
    const ZERO: Self;
    const ONE: Self;
    /// Bytes per element in the canonical little-endian encoding used for
    /// Merkle leaves and transcript absorption.
    const NUM_BYTES: usize;

    /// Reduce a `u64` into the field (used for AIR constants).
    fn from_u64(v: u64) -> Self;

    fn is_zero(self) -> bool {
        self == Self::ZERO
    }

    /// Multiplicative inverse. Panics on zero.
    fn inverse(self) -> Self;

    fn square(self) -> Self {
        self * self
    }

    fn double(self) -> Self {
        self + self
    }

    /// Exponentiation by squaring (for `u64` exponents).
    fn pow(self, mut exp: u64) -> Self {
        let mut base = self;
        let mut acc = Self::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                acc *= base;
            }
            base = base.square();
            exp >>= 1;
        }
        acc
    }

    /// Canonical little-endian encoding, `NUM_BYTES` long.
    fn to_bytes(self) -> Vec<u8>;

    /// Draw a uniform field element (Fiat–Shamir squeeze / masking).
    fn from_bytes_reduce(bytes: &[u8]) -> Self;

    fn random(rng: &mut impl Rng) -> Self;
}

/// An extension field `E` over base `B`: carries the base embedding and a
/// cheap base-scalar multiply, and can evaluate base/extension polynomials at
/// an extension point (Horner).
pub trait ExtField<B: Field>: Field + From<B> {
    /// Multiply by a base-field scalar (cheaper than a full extension multiply).
    fn mul_base(self, b: B) -> Self;

    /// Horner evaluation of a base-coefficient polynomial at this point.
    fn eval_base_poly(self, coeffs: &[B]) -> Self {
        let mut acc = Self::ZERO;
        for &c in coeffs.iter().rev() {
            acc = acc * self + Self::from(c);
        }
        acc
    }

    /// Horner evaluation of an extension-coefficient polynomial at this point.
    fn eval_ext_poly(self, coeffs: &[Self]) -> Self {
        let mut acc = Self::ZERO;
        for &c in coeffs.iter().rev() {
            acc = acc * self + c;
        }
        acc
    }
}

/// A field in which an AIR whose constants live in base `B` can be evaluated:
/// either `B` itself or an extension of it. Blanket-implemented for every
/// field over itself; each engine additionally implements it for its
/// extension field.
pub trait ConstraintField<B: Field>: Field + From<B> {}

impl<B: Field> ConstraintField<B> for B {}

// ── Field impls: classical (BabyBear) tower ────────────────────────────────

impl Field for BabyBear {
    const ZERO: Self = BabyBear { value: 0 };
    const ONE: Self = BabyBear { value: 1 };
    const NUM_BYTES: usize = 8;

    #[inline]
    fn from_u64(v: u64) -> Self {
        BabyBear::new(v)
    }
    #[inline]
    fn is_zero(self) -> bool {
        self.value == 0
    }
    #[inline]
    fn inverse(self) -> Self {
        BabyBear::inverse(self)
    }
    #[inline]
    fn pow(self, exp: u64) -> Self {
        BabyBear::pow(self, exp)
    }
    #[inline]
    fn to_bytes(self) -> Vec<u8> {
        BabyBear::to_bytes(&self).to_vec()
    }
    #[inline]
    fn from_bytes_reduce(bytes: &[u8]) -> Self {
        BabyBear::from_bytes_mod_order(bytes)
    }
    #[inline]
    fn random(rng: &mut impl Rng) -> Self {
        BabyBear::random(rng)
    }
}

impl Field for BabyBearExt {
    const ZERO: Self = BabyBearExt {
        c: [BabyBear { value: 0 }; 4],
    };
    const ONE: Self = BabyBearExt {
        c: [
            BabyBear { value: 1 },
            BabyBear { value: 0 },
            BabyBear { value: 0 },
            BabyBear { value: 0 },
        ],
    };
    const NUM_BYTES: usize = 32;

    #[inline]
    fn from_u64(v: u64) -> Self {
        BabyBearExt::from(BabyBear::new(v))
    }
    #[inline]
    fn is_zero(self) -> bool {
        BabyBearExt::is_zero(&self)
    }
    #[inline]
    fn inverse(self) -> Self {
        BabyBearExt::inverse(self)
    }
    #[inline]
    fn to_bytes(self) -> Vec<u8> {
        BabyBearExt::to_bytes(&self).to_vec()
    }
    #[inline]
    fn from_bytes_reduce(bytes: &[u8]) -> Self {
        // Four independent base squeezes.
        let mut c = [BabyBear::ZERO; 4];
        for (i, limb) in c.iter_mut().enumerate() {
            let start = (i * 8) % bytes.len().max(1);
            *limb = BabyBear::from_bytes_mod_order(&bytes[start..]);
        }
        BabyBearExt::new(c)
    }
    #[inline]
    fn random(rng: &mut impl Rng) -> Self {
        BabyBearExt::random(rng)
    }
}

impl ExtField<BabyBear> for BabyBearExt {
    #[inline]
    fn mul_base(self, b: BabyBear) -> Self {
        BabyBearExt::mul_base(self, b)
    }
}

impl ConstraintField<BabyBear> for BabyBearExt {}

// ── Field impls: circle (M31) tower ────────────────────────────────────────

#[cfg(feature = "circle")]
impl Field for M31 {
    const ZERO: Self = M31::ZERO;
    const ONE: Self = M31::ONE;
    const NUM_BYTES: usize = 4;

    #[inline]
    fn from_u64(v: u64) -> Self {
        M31::from_u64(v)
    }
    #[inline]
    fn is_zero(self) -> bool {
        M31::is_zero(self)
    }
    #[inline]
    fn inverse(self) -> Self {
        M31::inverse(self)
    }
    #[inline]
    fn pow(self, exp: u64) -> Self {
        M31::pow(self, exp)
    }
    #[inline]
    fn to_bytes(self) -> Vec<u8> {
        M31::to_bytes(self).to_vec()
    }
    #[inline]
    fn from_bytes_reduce(bytes: &[u8]) -> Self {
        M31::from_bytes_mod_order(bytes)
    }
    #[inline]
    fn random(rng: &mut impl Rng) -> Self {
        M31::random(rng)
    }
}

#[cfg(feature = "circle")]
impl Field for QM31 {
    const ZERO: Self = QM31::ZERO;
    const ONE: Self = QM31::ONE;
    const NUM_BYTES: usize = 16;

    #[inline]
    fn from_u64(v: u64) -> Self {
        QM31::from_base(M31::from_u64(v))
    }
    #[inline]
    fn is_zero(self) -> bool {
        QM31::is_zero(self)
    }
    #[inline]
    fn inverse(self) -> Self {
        QM31::inverse(self)
    }
    #[inline]
    fn to_bytes(self) -> Vec<u8> {
        QM31::to_bytes(self).to_vec()
    }
    #[inline]
    fn from_bytes_reduce(bytes: &[u8]) -> Self {
        let mut b = [0u8; 16];
        let n = bytes.len().min(16);
        b[..n].copy_from_slice(&bytes[..n]);
        QM31::from_bytes(&b)
    }
    #[inline]
    fn random(rng: &mut impl Rng) -> Self {
        QM31::random(rng)
    }
}

#[cfg(feature = "circle")]
impl ExtField<M31> for QM31 {
    #[inline]
    fn mul_base(self, b: M31) -> Self {
        QM31::mul_base(self, b)
    }
}

#[cfg(feature = "circle")]
impl ConstraintField<M31> for QM31 {}

// ── Montgomery batch inversion (generic) ───────────────────────────────────

/// Invert a slice with a single field inversion and `3(n−1)` multiplications.
/// Panics if any element is zero.
pub fn batch_inverse<F: Field>(values: &[F]) -> Vec<F> {
    if values.is_empty() {
        return Vec::new();
    }
    let mut prefix = Vec::with_capacity(values.len());
    let mut acc = F::ONE;
    for &v in values {
        prefix.push(acc);
        acc *= v;
    }
    let mut inv = acc.inverse();
    let mut out = vec![F::ZERO; values.len()];
    for i in (0..values.len()).rev() {
        out[i] = prefix[i] * inv;
        inv *= values[i];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn field_axioms<F: Field>(seed: u64) {
        let mut r = rand::rngs::StdRng::seed_from_u64(seed);
        for _ in 0..500 {
            let a = F::random(&mut r);
            let b = F::random(&mut r);
            let c = F::random(&mut r);
            assert_eq!(a + b, b + a);
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a - b + b, a);
            assert_eq!(a.double(), a + a);
            assert_eq!(a.square(), a * a);
            if !a.is_zero() {
                assert_eq!(a * a.inverse(), F::ONE);
            }
        }
    }

    #[test]
    fn babybear_is_a_field() {
        field_axioms::<BabyBear>(1);
        field_axioms::<BabyBearExt>(2);
    }

    #[cfg(feature = "circle")]
    #[test]
    fn m31_is_a_field() {
        field_axioms::<M31>(3);
        field_axioms::<QM31>(4);
    }

    #[test]
    fn batch_inverse_matches() {
        let vals: Vec<BabyBear> = (1..40u64).map(BabyBear::new).collect();
        for (v, inv) in vals.iter().zip(batch_inverse(&vals)) {
            assert_eq!(*v * inv, BabyBear::ONE);
        }
    }
}
