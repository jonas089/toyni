//! Composition polynomial for STARK proofs.
//!
//! Combines trace and constraints into H(x) = T(x) + Z_H(x) * sum(C_i(x)).

use ark_bls12_381::Fr;
use ark_ff::{Field, One, Zero};
use ark_poly::{
    EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial, univariate::DensePolynomial,
};

use crate::program::{constraints::ConstraintSystem, trace::ExecutionTrace};

/// Polynomial combining trace and constraints for STARK proofs.
pub struct CompositionPolynomial {
    /// The composed polynomial H(x)
    polynomial: DensePolynomial<Fr>,
    /// The evaluation domain
    domain: GeneralEvaluationDomain<Fr>,
}

impl CompositionPolynomial {
    /// Creates composition polynomial from trace and constraints.
    ///
    /// # Arguments
    ///
    /// * `trace` - The execution trace containing program state
    /// * `constraints` - The constraint system defining program rules
    /// * `domain` - The evaluation domain (can be extended)
    ///
    /// # Returns
    ///
    /// A new composition polynomial that encodes all constraints
    ///
    /// # Panics
    ///
    /// Panics if:
    /// * The domain size is not a power of 2
    /// * The trace height is greater than the domain size
    ///
    /// # Security Note
    ///
    /// The current implementation leaks information about the trace because:
    /// - The trace polynomial is directly included
    /// - No random masks are applied
    /// - The composition is not blinded
    pub fn new(
        trace: &ExecutionTrace,
        constraints: &ConstraintSystem,
        domain: GeneralEvaluationDomain<Fr>,
    ) -> Self {
        // Create constraint polynomial evaluations
        let mut constraint_evals = vec![Fr::zero(); domain.size()];

        // First, evaluate constraints on the original domain points
        let original_size = trace.height as usize;
        for (i, eval) in constraint_evals.iter_mut().enumerate().take(original_size) {
            let current_row = trace.get_column(i as u64);
            let next_row = trace.get_column(((i + 1) % original_size) as u64);

            for constraint in &constraints.transition_constraints {
                let constraint_eval = (constraint.evaluate)(current_row, next_row);
                *eval += constraint_eval;
            }
        }

        // Evaluate boundary constraints
        for constraint in &constraints.boundary_constraints {
            let row = trace.get_column(constraint.row);
            let eval = (constraint.evaluate)(row);
            constraint_evals[constraint.row as usize] += eval;
        }

        // If we have an extended domain, interpolate the constraint evaluations
        if domain.size() > original_size {
            // Create a polynomial from the original evaluations
            let original_domain = GeneralEvaluationDomain::new(original_size).unwrap();
            let constraint_poly = Evaluations::from_vec_and_domain(
                constraint_evals[..original_size].to_vec(),
                original_domain,
            )
            .interpolate();

            // Evaluate this polynomial over the extended domain
            constraint_evals = domain.fft(&constraint_poly.coeffs);
        }

        // Interpolate constraint polynomial sum(C_i(x))
        let constraint_poly =
            Evaluations::from_vec_and_domain(constraint_evals, domain).interpolate();

        // Create the vanishing polynomial Z_H(x)
        let mut z_h_evals = vec![Fr::one(); domain.size()];
        let omega = domain.element(1); // ω is the generator
        for (i, eval) in z_h_evals.iter_mut().enumerate().take(domain.size()) {
            let x = domain.element(i);
            for j in 0..domain.size() {
                let omega_j = omega.pow([j as u64]);
                *eval *= x - omega_j;
            }
        }
        let z_h = Evaluations::from_vec_and_domain(z_h_evals, domain).interpolate();

        // Compute H(x) = Z_H(x) * sum(C_i(x))
        let composition = &z_h * &constraint_poly;

        Self {
            polynomial: composition,
            domain,
        }
    }

    /// Creates polynomial from pre-computed evaluations.
    ///
    /// # Arguments
    ///
    /// * `evals` - The polynomial evaluations over the domain
    /// * `domain` - The evaluation domain
    ///
    /// # Returns
    ///
    /// A new composition polynomial with the given evaluations
    pub fn from_evaluations(evals: Vec<Fr>, domain: GeneralEvaluationDomain<Fr>) -> Self {
        let poly = Evaluations::from_vec_and_domain(evals.clone(), domain).interpolate();
        Self {
            polynomial: poly,
            domain,
        }
    }

    /// Returns polynomial degree.
    ///
    /// The degree is important for:
    /// - FRI low-degree testing
    /// - Proof size estimation
    /// - Security parameter selection
    ///
    /// # Returns
    ///
    /// The degree of the polynomial
    pub fn degree(&self) -> usize {
        self.polynomial.degree()
    }

    /// Evaluates polynomial at point.
    ///
    /// # Arguments
    ///
    /// * `point` - The point at which to evaluate
    ///
    /// # Returns
    ///
    /// The value of the polynomial at the given point
    pub fn evaluate(&self, point: Fr) -> Fr {
        self.polynomial.evaluate(&point)
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
        &self.polynomial.coeffs
    }

    /// Returns evaluations over domain.
    ///
    /// # Returns
    ///
    /// A vector of evaluations at each point in the domain
    ///
    /// # Security Note
    ///
    /// The current implementation exposes all evaluations, which may leak
    /// information about the trace. In a zero-knowledge implementation, these
    /// should be committed to using a Merkle tree.
    pub fn evaluations(&self) -> Vec<Fr> {
        self.domain.fft(&self.polynomial.coeffs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_composition_polynomial() {
        // Create a simple trace: x[n] = n
        let mut trace = ExecutionTrace::new(4, 1);
        for i in 0..4 {
            let mut column = HashMap::new();
            column.insert("x".to_string(), i);
            trace.insert_column(column);
        }

        // Create constraint system: x[n] = x[n-1] + 1
        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_current = current.get("x").unwrap();
                let x_next = next.get("x").unwrap();
                Fr::from(*x_next) - Fr::from(*x_current + 1)
            }),
        );

        // Create evaluation domain
        let domain = GeneralEvaluationDomain::<Fr>::new(4).unwrap();

        // Create composition polynomial
        let comp_poly = CompositionPolynomial::new(&trace, &constraints, domain);

        // The composition polynomial should evaluate to zero at all points where:
        // 1. The trace values are correct (x[n] = n)
        // 2. The constraints are satisfied (x[n] = x[n-1] + 1)
        let evals = comp_poly.evaluations();
        for eval in evals {
            assert!(eval.is_zero());
        }
    }
}
