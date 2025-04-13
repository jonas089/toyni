//! Basic polynomial operations over finite fields.

use ark_bls12_381::Fr;
use ark_ff::Zero;
use ark_poly::univariate::DensePolynomial;
use std::fmt;

#[derive(Debug, Clone)]
/// Polynomial with finite field coefficients.
///
/// The polynomial is stored as a vector of coefficients, where the index represents
/// the power of x. For example, [1, 2, 3] represents 3x² + 2x + 1.
///
/// # Invariants
///
/// * The coefficients vector should not have trailing zeros
/// * All coefficients should be valid field elements
pub struct Polynomial {
    /// Coefficients in ascending order of power.
    /// The vector must not have trailing zeros.
    pub coefficients: Vec<Fr>,
}

impl Polynomial {
    /// Creates a new polynomial from coefficients.
    ///
    /// # Arguments
    ///
    /// * `coefficients` - The coefficients of the polynomial
    ///
    /// # Returns
    ///
    /// A new polynomial with the given coefficients
    ///
    /// # Note
    ///
    /// Any trailing zeros in the coefficients vector will be removed.
    pub fn new(mut coefficients: Vec<Fr>) -> Self {
        // Remove trailing zeros
        while coefficients.last().is_some_and(|&x| x.is_zero()) {
            coefficients.pop();
        }
        Self { coefficients }
    }

    /// Returns the degree of the polynomial.
    ///
    /// The degree is the highest power of x with a non-zero coefficient.
    /// For the zero polynomial, the degree is 0.
    ///
    /// # Returns
    ///
    /// The degree of the polynomial
    pub fn degree(&self) -> usize {
        if self.coefficients.is_empty() {
            0
        } else {
            self.coefficients.len() - 1
        }
    }

    /// Returns the coefficient of the highest power term.
    ///
    /// The leading coefficient is the coefficient of the highest power term.
    /// For the zero polynomial, returns 0.
    ///
    /// # Returns
    ///
    /// The leading coefficient
    pub fn leading_coefficient(&self) -> Fr {
        self.coefficients.last().copied().unwrap_or(Fr::zero())
    }

    /// Divides this polynomial by another, returns (quotient, remainder).
    ///
    /// # Arguments
    ///
    /// * `divisor` - The polynomial to divide by
    ///
    /// # Returns
    ///
    /// Option containing a tuple of (quotient, remainder) if division is possible,
    /// None if the divisor is zero or empty
    ///
    /// # Details
    ///
    /// The division is performed using the standard long division algorithm
    /// over the finite field. The remainder will have degree less than the divisor.
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

    /// Adds two polynomials.
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

    /// Multiplies two polynomials.
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

    /// Evaluates the polynomial at point x.
    ///
    /// # Arguments
    ///
    /// * `x` - The point at which to evaluate the polynomial
    ///
    /// # Returns
    ///
    /// The value of the polynomial at x
    ///
    /// # Details
    ///
    /// Evaluation is performed using Horner's method for efficiency.
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

    /// Creates a polynomial from a dense polynomial.
    ///
    /// # Arguments
    ///
    /// * `poly` - The dense polynomial to convert
    ///
    /// # Returns
    ///
    /// A new polynomial with the same coefficients
    pub fn from_dense_poly(poly: DensePolynomial<Fr>) -> Self {
        Self::new(poly.coeffs)
    }

    /// Checks if the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|c| c.is_zero())
    }

    /// Creates the zero polynomial.
    ///
    /// # Returns
    ///
    /// A polynomial representing zero
    pub fn zero() -> Self {
        Self::new(vec![Fr::zero()])
    }

    /// Returns polynomial coefficients.
    ///
    /// # Returns
    ///
    /// A slice containing the polynomial coefficients
    ///
    /// # Security Note
    ///
    /// The current implementation exposes the raw coefficients, which may leak
    /// information about the trace. In a zero-knowledge implementation, these
    /// should be committed to using a Merkle tree.
    pub fn coefficients(&self) -> &[Fr] {
        &self.coefficients
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
