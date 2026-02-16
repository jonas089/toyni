//! Basic polynomial operations over BabyBear field.

use crate::babybear::BabyBear;
use std::fmt;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coefficients: Vec<BabyBear>,
}

impl Polynomial {
    pub fn new(mut coefficients: Vec<BabyBear>) -> Self {
        while coefficients.last().is_some_and(|x| x.is_zero()) {
            coefficients.pop();
        }
        Self { coefficients }
    }

    pub fn degree(&self) -> usize {
        if self.coefficients.is_empty() {
            0
        } else {
            self.coefficients.len() - 1
        }
    }

    pub fn leading_coefficient(&self) -> BabyBear {
        self.coefficients
            .last()
            .copied()
            .unwrap_or(BabyBear::zero())
    }

    pub fn divide(&self, divisor: &Polynomial) -> Option<(Polynomial, Polynomial)> {
        if divisor.coefficients.is_empty() || divisor.leading_coefficient().is_zero() {
            return None;
        }

        let dividend = self.coefficients.clone();
        let divisor_degree = divisor.degree();
        let dividend_degree = self.degree();

        if dividend_degree < divisor_degree {
            return Some((Polynomial::zero(), self.clone()));
        }

        let mut quotient = vec![BabyBear::zero(); dividend_degree - divisor_degree + 1];
        let mut remainder = dividend.clone();

        for i in (0..=dividend_degree - divisor_degree).rev() {
            let leading_coeff = remainder[i + divisor_degree];
            if leading_coeff.is_zero() {
                continue;
            }

            quotient[i] = leading_coeff / divisor.leading_coefficient();

            for j in 0..=divisor_degree {
                remainder[i + j] = remainder[i + j] - quotient[i] * divisor.coefficients[j];
            }
        }

        while !remainder.is_empty() && remainder.last().unwrap().is_zero() {
            remainder.pop();
        }

        Some((Polynomial::new(quotient), Polynomial::new(remainder)))
    }

    pub fn divide_by_linear(&self, z: BabyBear) -> (Polynomial, BabyBear) {
        let mut quotient = vec![BabyBear::zero(); self.coefficients.len().saturating_sub(1)];
        let mut remainder = BabyBear::zero();

        let mut acc = BabyBear::zero();
        for (i, &coeff) in self.coefficients.iter().rev().enumerate() {
            if i == 0 {
                remainder = acc;
                break;
            }
            acc = coeff + z * acc;
            let length = quotient.len();
            quotient[length - i] = acc;
        }

        (Polynomial::new(quotient), remainder)
    }

    pub fn add(&self, other: &Polynomial) -> Polynomial {
        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut result = vec![BabyBear::zero(); max_len];

        for (i, coeff) in result.iter_mut().enumerate().take(self.coefficients.len()) {
            *coeff = *coeff + self.coefficients[i];
        }

        for (i, coeff) in result.iter_mut().enumerate().take(other.coefficients.len()) {
            *coeff = *coeff + other.coefficients[i];
        }

        Polynomial::new(result)
    }

    pub fn sub(&self, other: &Polynomial) -> Polynomial {
        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut result = vec![BabyBear::zero(); max_len];

        for i in 0..self.coefficients.len() {
            result[i] = self.coefficients[i];
        }

        for i in 0..other.coefficients.len() {
            result[i] = result[i] - other.coefficients[i];
        }

        Polynomial::new(result)
    }

    pub fn multiply(&self, other: &Polynomial) -> Polynomial {
        if self.coefficients.is_empty() || other.coefficients.is_empty() {
            return Polynomial::new(vec![]);
        }

        let mut result = vec![BabyBear::zero(); self.degree() + other.degree() + 1];

        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                result[i + j] = result[i + j] + self.coefficients[i] * other.coefficients[j];
            }
        }

        Polynomial::new(result)
    }

    pub fn evaluate(&self, x: BabyBear) -> BabyBear {
        if self.coefficients.is_empty() {
            return BabyBear::zero();
        }

        let mut result = self.coefficients[self.coefficients.len() - 1];
        for &coeff in self.coefficients.iter().rev().skip(1) {
            result = result * x + coeff;
        }
        result
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|c| c.is_zero())
    }

    pub fn zero() -> Self {
        Self::new(vec![BabyBear::zero()])
    }

    pub fn mul(&self, other: &Self) -> Self {
        self.multiply(other)
    }

    pub fn coefficients(&self) -> &[BabyBear] {
        &self.coefficients
    }

    pub fn scale(&self, scalar: BabyBear) -> Self {
        let scaled_coeffs: Vec<BabyBear> = self.coefficients.iter().map(|c| *c * scalar).collect();
        Polynomial::new(scaled_coeffs)
    }

    pub fn lagrange_interpolate(xs: &[BabyBear], ys: &[BabyBear]) -> Self {
        assert_eq!(xs.len(), ys.len(), "Mismatched input lengths");
        let n = xs.len();
        let mut result = Polynomial::zero();

        for i in 0..n {
            let mut basis = Polynomial::new(vec![BabyBear::one()]);
            let mut denom = BabyBear::one();

            for j in 0..n {
                if i == j {
                    continue;
                }
                let factor = Polynomial::new(vec![-xs[j], BabyBear::one()]);
                basis = basis.mul(&factor);
                denom = denom * (xs[i] - xs[j]);
            }

            let coeff = ys[i] / denom;
            basis = basis.scale(coeff);
            result = result.add(&basis);
        }

        result
    }
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.coefficients.is_empty() {
            return write!(f, "0");
        }

        let mut terms = Vec::new();
        for (i, &coeff) in self.coefficients.iter().enumerate() {
            if !coeff.is_zero() {
                let term = if i == 0 {
                    format!("{}", coeff)
                } else if i == 1 {
                    format!("{}x", coeff)
                } else {
                    format!("{}x^{}", coeff, i)
                };
                terms.push(term);
            }
        }

        if terms.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", terms.join(" + "))
        }
    }
}
