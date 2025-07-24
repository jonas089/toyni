//! Basic polynomial operations over finite fields.

use ark_bls12_381::Fr;
use ark_ff::{UniformRand, Zero};
use ark_poly::univariate::DensePolynomial;
use rand;
use std::fmt;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coefficients: Vec<Fr>,
}

impl Polynomial {
    pub fn new(mut coefficients: Vec<Fr>) -> Self {
        while coefficients.last().is_some_and(|&x| x.is_zero()) {
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

    pub fn leading_coefficient(&self) -> Fr {
        self.coefficients.last().copied().unwrap_or(Fr::zero())
    }

    pub fn divide(&self, divisor: &Polynomial) -> Option<(Polynomial, Polynomial)> {
        if divisor.coefficients.is_empty() || divisor.leading_coefficient().is_zero() {
            return None;
        }

        let dividend = self.coefficients.clone();
        let divisor_degree = divisor.degree();
        let dividend_degree = self.degree();

        // If dividend degree is less than divisor degree, quotient is zero
        if dividend_degree < divisor_degree {
            return Some((Polynomial::zero(), self.clone()));
        }

        let mut quotient = vec![Fr::zero(); dividend_degree - divisor_degree + 1];
        let mut remainder = dividend.clone();

        // Perform long division
        for i in (0..=dividend_degree - divisor_degree).rev() {
            let leading_coeff = remainder[i + divisor_degree];
            if leading_coeff.is_zero() {
                continue;
            }

            quotient[i] = leading_coeff / divisor.leading_coefficient();

            // Subtract divisor * quotient term from remainder
            for j in 0..=divisor_degree {
                remainder[i + j] -= quotient[i] * divisor.coefficients[j];
            }
        }

        // Trim leading zeros from remainder
        while !remainder.is_empty() && remainder.last().unwrap().is_zero() {
            remainder.pop();
        }

        Some((Polynomial::new(quotient), Polynomial::new(remainder)))
    }

    pub fn add(&self, other: &Polynomial) -> Polynomial {
        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut result = vec![Fr::zero(); max_len];

        for (i, coeff) in result.iter_mut().enumerate().take(self.coefficients.len()) {
            *coeff += self.coefficients[i];
        }

        for (i, coeff) in result.iter_mut().enumerate().take(other.coefficients.len()) {
            *coeff += other.coefficients[i];
        }

        Polynomial::new(result)
    }

    pub fn sub(&self, other: &Polynomial) -> Polynomial {
        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut result = vec![Fr::zero(); max_len];

        for i in 0..self.coefficients.len() {
            result[i] = self.coefficients[i];
        }

        for i in 0..other.coefficients.len() {
            result[i] -= other.coefficients[i];
        }

        Polynomial::new(result)
    }

    pub fn multiply(&self, other: &Polynomial) -> Polynomial {
        if self.coefficients.is_empty() || other.coefficients.is_empty() {
            return Polynomial::new(vec![]);
        }

        let mut result = vec![Fr::zero(); self.degree() + other.degree() + 1];

        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                result[i + j] += self.coefficients[i] * other.coefficients[j];
            }
        }

        Polynomial::new(result)
    }

    pub fn evaluate(&self, x: Fr) -> Fr {
        if self.coefficients.is_empty() {
            return Fr::zero();
        }

        let mut result = self.coefficients[self.coefficients.len() - 1];
        for &coeff in self.coefficients.iter().rev().skip(1) {
            result = result * x + coeff;
        }
        result
    }

    pub fn evaluate_at_u64(&self, x: u64) -> Fr {
        self.evaluate(Fr::from(x))
    }

    pub fn from_dense_poly(poly: DensePolynomial<Fr>) -> Self {
        Self::new(poly.coeffs)
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|c| c.is_zero())
    }

    pub fn zero() -> Self {
        Self::new(vec![Fr::zero()])
    }

    pub fn random(degree: usize, rng: &mut impl rand::Rng) -> Self {
        let mut coefficients = Vec::with_capacity(degree + 1);
        for _ in 0..=degree {
            coefficients.push(Fr::rand(rng));
        }
        Self::new(coefficients)
    }

    pub fn mul(&self, other: &Self) -> Self {
        self.multiply(other)
    }

    pub fn coefficients(&self) -> &[Fr] {
        &self.coefficients
    }

    pub fn scale(&self, scalar: Fr) -> Self {
        let scaled_coeffs: Vec<Fr> = self.coefficients.iter().map(|c| *c * scalar).collect();
        Polynomial::new(scaled_coeffs)
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
