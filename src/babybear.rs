// BabyBear field implementation
// Prime: p = 2^31 - 2^27 + 1 = 2013265921
// Ported from joda for use in toyni's STARK prover

use std::hash::{Hash, Hasher};
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

pub const BABYBEAR_PRIME: u32 = 2013265921; // 2^31 - 2^27 + 1

/// A BabyBear field element, always kept in canonical form (`0 <= value < p`),
/// which fits in a `u32`. `repr(C)` with a single `u32` field means a
/// `&[BabyBear]` has the same layout as a `&[u32]` and can be handed to the GPU
/// NTT without conversion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
pub struct BabyBear {
    pub value: u32,
}

impl Hash for BabyBear {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.value.hash(state);
    }
}

impl BabyBear {
    pub const PRIME: u32 = BABYBEAR_PRIME;

    #[inline]
    pub fn new(value: u64) -> Self {
        Self {
            value: (value % Self::PRIME as u64) as u32,
        }
    }

    #[inline]
    pub fn zero() -> Self {
        Self { value: 0 }
    }

    #[inline]
    pub fn one() -> Self {
        Self { value: 1 }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.value == 0
    }

    #[inline]
    pub fn from_u32(value: u32) -> Self {
        Self::new(value as u64)
    }

    #[inline]
    pub fn to_bytes(&self) -> [u8; 4] {
        self.value.to_le_bytes()
    }

    #[inline]
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut arr = [0u8; 4];
        arr.copy_from_slice(&bytes[..4]);
        Self::new(u32::from_le_bytes(arr) as u64)
    }

    /// Create a BabyBear element from a byte slice, reducing mod p.
    ///
    /// Still folds 8 bytes even though an element is only 4: reducing a 64-bit
    /// draw leaves negligible modular bias, whereas a 32-bit draw would be
    /// biased by roughly 2:1 across residues (2^32 / p ~ 2.13).
    pub fn from_bytes_mod_order(bytes: &[u8]) -> Self {
        let mut val: u64 = 0;
        for (i, &byte) in bytes.iter().take(8).enumerate() {
            val |= (byte as u64) << (i * 8);
        }
        Self::new(val)
    }

    /// Generate a random BabyBear element.
    pub fn random(rng: &mut impl rand::Rng) -> Self {
        let val: u32 = rng.gen_range(0..Self::PRIME);
        Self { value: val }
    }

    /// Fold a value in [0, 2p) into [0, p). Both operands of an add are
    /// canonical, so their sum is < 2p < 2^32 and never overflows a u32.
    #[inline]
    fn reduce_sum(val: u32) -> u32 {
        if val >= Self::PRIME { val - Self::PRIME } else { val }
    }

    pub fn pow(self, mut exp: u64) -> Self {
        if exp == 0 {
            return Self::one();
        }

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

    /// Modular inverse using Fermat's little theorem: a^(p-2) mod p
    pub fn inverse(self) -> Self {
        assert_ne!(self.value, 0, "Cannot invert zero");
        self.pow(Self::PRIME as u64 - 2)
    }

    /// Get the primitive root of unity for a given power of 2.
    /// BabyBear supports NTT up to size 2^27.
    pub fn get_root_of_unity(log_n: u32) -> Self {
        assert!(log_n <= 27, "BabyBear only supports NTT up to 2^27");

        const TWO_ADICITY: u32 = 27;
        const PRIMITIVE_ROOT_OF_UNITY: u64 = 440564289; // 31^15 mod p

        let exp = 1u64 << (TWO_ADICITY - log_n);
        Self::new(PRIMITIVE_ROOT_OF_UNITY).pow(exp)
    }
}

impl Add for BabyBear {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        let sum = self.value + rhs.value;
        Self {
            value: Self::reduce_sum(sum),
        }
    }
}

impl AddAssign for BabyBear {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for BabyBear {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let diff = if self.value >= rhs.value {
            self.value - rhs.value
        } else {
            self.value + Self::PRIME - rhs.value
        };
        Self { value: diff }
    }
}

impl SubAssign for BabyBear {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for BabyBear {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let product = (self.value as u64) * (rhs.value as u64);
        let reduced = (product % Self::PRIME as u64) as u32;
        Self { value: reduced }
    }
}

impl MulAssign for BabyBear {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Div for BabyBear {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl Neg for BabyBear {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        if self.value == 0 {
            self
        } else {
            Self {
                value: Self::PRIME - self.value,
            }
        }
    }
}

impl std::fmt::Display for BabyBear {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_arithmetic() {
        let a = BabyBear::new(100);
        let b = BabyBear::new(200);

        let sum = a + b;
        assert_eq!(sum.value, 300);

        let diff = b - a;
        assert_eq!(diff.value, 100);

        let prod = a * b;
        assert_eq!(prod.value, 20000);
    }

    #[test]
    fn test_modular_reduction() {
        let a = BabyBear::new(BABYBEAR_PRIME as u64 + 5);
        assert_eq!(a.value, 5);
    }

    #[test]
    fn test_inverse() {
        let a = BabyBear::new(7);
        let inv = a.inverse();
        let prod = a * inv;
        assert_eq!(prod.value, 1);
    }

    #[test]
    fn test_pow() {
        let a = BabyBear::new(3);
        let result = a.pow(4);
        assert_eq!(result.value, 81);
    }

    #[test]
    fn test_root_of_unity() {
        for log_n in 1..=10 {
            let omega = BabyBear::get_root_of_unity(log_n);
            let n = 1u64 << log_n;
            let result = omega.pow(n);
            assert_eq!(
                result.value, 1,
                "omega^n should equal 1 for log_n={}",
                log_n
            );
        }
    }

    #[test]
    fn test_negation() {
        let a = BabyBear::new(100);
        let neg_a = -a;
        let sum = a + neg_a;
        assert_eq!(sum.value, 0);
    }

    /// The u32 representation is only sound if no intermediate overflows.
    /// Canonical operands are < p < 2^31, so an add reaches at most 2p-2 and
    /// `Sub`/`Neg`'s `value + PRIME` at most 2p-1 — both under 2^32. These run
    /// at the boundary, and a debug build panics on overflow.
    #[test]
    fn test_no_overflow_at_boundary() {
        let max = BabyBear::new(BABYBEAR_PRIME as u64 - 1); // p - 1, the largest element
        assert_eq!((max + max).value, BABYBEAR_PRIME - 2);  // sum 2p-2 folds back
        assert_eq!((BabyBear::zero() - max).value, 1);      // borrow path
        assert_eq!((-max).value, 1);
        assert_eq!((max * max).value, 1);                   // (-1)^2 = 1
        assert_eq!(max.pow(2).value, 1);
        assert_eq!((max + BabyBear::one()).value, 0);
    }

    /// Every operation must leave the value canonical, since the GPU NTT and the
    /// byte encoding both assume it.
    #[test]
    fn test_results_stay_canonical() {
        let mut r = {
            use rand::SeedableRng;
            rand::rngs::StdRng::seed_from_u64(7)
        };
        for _ in 0..2000 {
            let a = BabyBear::random(&mut r);
            let b = BabyBear::random(&mut r);
            for v in [a + b, a - b, a * b, -a, a.pow(12345)] {
                assert!(v.value < BABYBEAR_PRIME, "non-canonical: {}", v.value);
            }
        }
    }

    /// A slice of BabyBear must be layout-compatible with a slice of u32; the
    /// CUDA path reinterprets one as the other with no copy.
    #[test]
    fn test_u32_layout() {
        assert_eq!(std::mem::size_of::<BabyBear>(), std::mem::size_of::<u32>());
        assert_eq!(std::mem::align_of::<BabyBear>(), std::mem::align_of::<u32>());
    }

    #[test]
    fn test_byte_roundtrip() {
        for v in [0u64, 1, 12345, BABYBEAR_PRIME as u64 - 1] {
            let a = BabyBear::new(v);
            assert_eq!(a.to_bytes().len(), 4);
            assert_eq!(BabyBear::from_bytes(&a.to_bytes()), a);
        }
    }

    #[test]
    fn test_division() {
        let a = BabyBear::new(100);
        let b = BabyBear::new(7);
        let q = a / b;
        assert_eq!((q * b).value, a.value);
    }
}
