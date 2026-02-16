// BabyBear field implementation
// Prime: p = 2^31 - 2^27 + 1 = 2013265921
// Ported from joda for use in toyni's STARK prover

use std::hash::{Hash, Hasher};
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

pub const BABYBEAR_PRIME: u64 = 2013265921; // 2^31 - 2^27 + 1

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
pub struct BabyBear {
    pub value: u64,
}

impl Hash for BabyBear {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.value.hash(state);
    }
}

impl BabyBear {
    pub const PRIME: u64 = BABYBEAR_PRIME;

    #[inline]
    pub fn new(value: u64) -> Self {
        Self {
            value: value % Self::PRIME,
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
    pub fn to_bytes(&self) -> [u8; 8] {
        self.value.to_le_bytes()
    }

    #[inline]
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut arr = [0u8; 8];
        arr.copy_from_slice(&bytes[..8]);
        Self::new(u64::from_le_bytes(arr))
    }

    /// Create a BabyBear element from a byte slice, reducing mod p.
    pub fn from_bytes_mod_order(bytes: &[u8]) -> Self {
        let mut val: u64 = 0;
        for (i, &byte) in bytes.iter().take(8).enumerate() {
            val |= (byte as u64) << (i * 8);
        }
        Self::new(val)
    }

    /// Generate a random BabyBear element.
    pub fn random(rng: &mut impl rand::Rng) -> Self {
        let val: u64 = rng.gen_range(0..Self::PRIME);
        Self { value: val }
    }

    #[inline]
    fn reduce_wide(val: u64) -> u64 {
        let mut v = val;
        if v >= Self::PRIME {
            v -= Self::PRIME;
        }
        if v >= Self::PRIME {
            v -= Self::PRIME;
        }
        v
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
        self.pow(Self::PRIME - 2)
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
            value: Self::reduce_wide(sum),
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
        let product = (self.value as u128) * (rhs.value as u128);
        let reduced = (product % Self::PRIME as u128) as u64;
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
        let a = BabyBear::new(BABYBEAR_PRIME + 5);
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

    #[test]
    fn test_division() {
        let a = BabyBear::new(100);
        let b = BabyBear::new(7);
        let q = a / b;
        assert_eq!((q * b).value, a.value);
    }
}
