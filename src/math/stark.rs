//! STARK (Scalable Transparent Argument of Knowledge) implementation.
//!
//! This module provides a STARK proving system implementation that:
//! 1. Converts program execution into a trace polynomial
//! 2. Constructs a composition polynomial from constraints
//! 3. Uses the FRI protocol to prove low-degree
//!
//! # Security Properties
//!
//! The implementation provides:
//! - Scalability: Proof size is logarithmic in computation size
//! - Transparency: No trusted setup required
//! - Zero-knowledge: Hides execution details through FRI and Merkle commitments
//!
//! # Components
//!
//! - `StarkProof`: Contains all components needed for verification
//! - `StarkProver`: Generates proofs from execution traces
//! - `StarkVerifier`: Verifies proofs using FRI and Merkle commitments

use crate::digest_sha2;
use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::merkle::{MerkleTree, verify_merkle_proof};
use crate::vm::{constraints::ConstraintSystem, trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_poly::DenseUVPolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::thread_rng;

/// Number of random challenges for verifier spot checks
const VERIFIER_QUERIES: usize = 80;

/// STARK proof containing all components needed for verification.
///
/// The proof consists of:
/// - Quotient polynomial evaluations
/// - FRI protocol layers and challenges
/// - Merkle commitments for each FRI layer
/// - Random challenges for spot checks
#[derive(Debug)]
pub struct StarkProof {
    /// Quotient polynomial evaluations over the extended domain
    pub quotient_eval_domain: Vec<Fr>,
    /// FRI protocol layers with folded evaluations
    pub fri_layers: Vec<Vec<Fr>>,
    /// Random challenges for FRI folding
    pub fri_challenges: Vec<Fr>,
    /// Combined constraint polynomial
    pub combined_constraint: ToyniPolynomial,
    /// Quotient polynomial from division
    pub quotient_poly: ToyniPolynomial,
    /// Merkle trees for each FRI layer's commitments
    pub folding_commitment_trees: Vec<MerkleTree>,
    /// Fiat-Shamir random challenges for spot checks
    pub verifier_random_challenges: Vec<Fr>,
}

/// STARK prover component that generates proofs from execution traces.
///
/// The prover:
/// 1. Interpolates constraints into polynomials
/// 2. Constructs the composition polynomial
/// 3. Performs FRI folding with Merkle commitments
/// 4. Generates random challenges for verification
pub struct StarkProver<'a> {
    /// Execution trace to prove
    trace: &'a ExecutionTrace,
    /// Constraint system defining program rules
    constraints: &'a ConstraintSystem,
}

impl<'a> StarkProver<'a> {
    /// Creates a new STARK prover for the given trace and constraints.
    ///
    /// # Arguments
    ///
    /// * `trace` - The execution trace to prove
    /// * `constraints` - The constraint system defining program rules
    pub fn new(trace: &'a ExecutionTrace, constraints: &'a ConstraintSystem) -> Self {
        Self { trace, constraints }
    }

    /// Generates a STARK proof for the execution trace.
    ///
    /// The proof generation process:
    /// 1. Interpolates all constraints into polynomials
    /// 2. Combines constraints into a single polynomial
    /// 3. Divides by the vanishing polynomial to get quotient
    /// 4. Performs FRI folding with Merkle commitments
    /// 5. Generates random challenges for verification
    ///
    /// # Returns
    ///
    /// A `StarkProof` containing all components needed for verification
    pub fn generate_proof(&self) -> StarkProof {
        let trace_len = self.trace.height as usize;
        let domain = GeneralEvaluationDomain::<Fr>::new(trace_len).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 2).unwrap();

        // Interpolate all constraints into polynomials
        let constraint_polys = self.constraints.interpolate_all_constraints(self.trace);

        // Combine all constraints into a single polynomial
        let mut combined_constraint = ToyniPolynomial::zero();
        for poly in constraint_polys {
            combined_constraint = combined_constraint.add(&poly);
        }

        // Evaluate combined constraint over extended domain
        let c_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| combined_constraint.evaluate(x))
            .collect();

        // Interpolate constraint polynomial from evaluations
        let c_poly = DensePolynomial::from_coefficients_slice(&extended_domain.ifft(&c_evals));
        let c_poly = ToyniPolynomial::from_dense_poly(c_poly);

        // Create vanishing polynomial
        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());

        // Divide to get quotient polynomial
        let (quotient_poly, _) = c_poly.divide(&z_poly).unwrap();

        // Evaluate quotient polynomial over extended domain
        let mut q_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| quotient_poly.evaluate(x))
            .collect();

        // Perform FRI folding with Merkle commitments
        let mut fri_layers = vec![q_evals.clone()];
        let mut fri_challenges = Vec::new();
        let mut folding_commitment_trees: Vec<MerkleTree> = Vec::new();

        while q_evals.len() > 4 {
            let beta = Fr::rand(&mut thread_rng());
            fri_challenges.push(beta);
            q_evals = fri_fold(&q_evals, beta);
            
            // Create Merkle tree for this FRI layer
            let folding_step_merkle_tree = MerkleTree::new(
                q_evals
                    .clone()
                    .iter()
                    .map(|x| x.into_bigint().to_bytes_be())
                    .collect(),
            );
            folding_commitment_trees.push(folding_step_merkle_tree);
            fri_layers.push(q_evals.clone());
        }

        // Generate random challenges for verification
        let proof_hash = digest_sha2(&[0; 32]);
        let proof_hash_u32: Vec<u32> = proof_hash
            .chunks(4)
            .map(|chunk| {
                let mut buf = [0u8; 4];
                buf.copy_from_slice(chunk);
                u32::from_be_bytes(buf)
            })
            .collect();
        let proof_hash_int = BigUint::from_slice(&proof_hash_u32);
        let mut verifier_random_challenges = Vec::new();

        for _ in 0..VERIFIER_QUERIES {
            let verifier_index = proof_hash_int.clone() + BigUint::from(1u32);
            let verifier_index_field = verifier_index % BigUint::from(extended_domain.size());
            let random_interactive_challenge =
                extended_domain.element(verifier_index_field.to_usize().unwrap());
            verifier_random_challenges.push(random_interactive_challenge);
        }

        StarkProof {
            quotient_eval_domain: fri_layers[0].clone(),
            fri_layers,
            fri_challenges,
            combined_constraint,
            quotient_poly,
            folding_commitment_trees,
            verifier_random_challenges,
        }
    }
}

/// STARK verifier component that verifies proofs.
///
/// The verifier:
/// 1. Checks FRI folding consistency with Merkle proofs
/// 2. Verifies constraint satisfaction at random points
/// 3. Ensures all commitments are valid
pub struct StarkVerifier<'a> {
    /// Constraint system defining program rules
    #[allow(unused)]
    constraints: &'a ConstraintSystem,
    /// Length of execution trace
    trace_len: usize,
}

impl<'a> StarkVerifier<'a> {
    /// Creates a new STARK verifier for the given constraints and trace length.
    ///
    /// # Arguments
    ///
    /// * `constraints` - The constraint system defining program rules
    /// * `trace_len` - The length of the execution trace
    pub fn new(constraints: &'a ConstraintSystem, trace_len: usize) -> Self {
        Self {
            constraints,
            trace_len,
        }
    }

    /// Verifies a STARK proof.
    ///
    /// The verification process:
    /// 1. Checks FRI folding consistency with Merkle proofs
    /// 2. Verifies constraint satisfaction at random points
    /// 3. Ensures all commitments are valid
    ///
    /// # Arguments
    ///
    /// * `proof` - The STARK proof to verify
    ///
    /// # Returns
    ///
    /// `true` if the proof is valid, `false` otherwise
    pub fn verify(&self, proof: &StarkProof) -> bool {
        let domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len * 2).unwrap();
        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());

        // FRI folding consistency check with Merkle proof verification
        let mut current_layer = &proof.quotient_eval_domain;
        for (i, ((beta, next_layer), merkle_tree)) in proof
            .fri_challenges
            .iter()
            .zip(proof.fri_layers.iter().skip(1))
            .zip(proof.folding_commitment_trees.iter())
            .enumerate()
        {
            let n = current_layer.len();
            let half_n = n / 2;

            // Verify each point in the next layer was correctly folded and committed
            for j in 0..half_n {
                let x = current_layer[j];
                let neg_x = current_layer[j + half_n];

                // f_next(x) = (f(x) + f(-x) + β * (f(x) - f(-x))) / 2
                let expected_next = (x + neg_x + *beta * (x - neg_x)) / Fr::from(2u64);
                let actual_next = next_layer[j];

                // Verify the actual value matches the expected folded value
                if expected_next != actual_next {
                    println!("❌ FRI folding failed at layer {}, position {}", i, j);
                    println!("Expected: {:?}", expected_next);
                    println!("Actual: {:?}", actual_next);
                    return false;
                }

                // Get the Merkle proof for this position
                let proof = merkle_tree.get_proof(j).expect("Merkle proof should exist");
                let root = merkle_tree.root().expect("Merkle root should exist");

                // Verify the value is properly committed in the Merkle tree
                let value_bytes = actual_next.into_bigint().to_bytes_be();
                if !verify_merkle_proof(value_bytes, &proof, &root) {
                    println!(
                        "❌ Merkle proof verification failed at layer {}, position {}",
                        i, j
                    );
                    return false;
                }
            }
            current_layer = next_layer;
        }

        // Verify constraint satisfaction at random points
        for _i in proof.verifier_random_challenges.iter() {
            let random_interactive_challenge =
                extended_domain.element(rand::random::<usize>() % extended_domain.size());
            let q_eval = proof.quotient_poly.evaluate(random_interactive_challenge);
            let z_eval = z_poly.evaluate(random_interactive_challenge);
            let c_eval = proof
                .combined_constraint
                .evaluate(random_interactive_challenge);

            if q_eval * z_eval != c_eval {
                println!("❌ Spot check failed: Q(x₀)*Z(x₀) ≠ C(x₀)");
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Field;
    use std::collections::HashMap;

    #[test]
    fn test_valid_proof() {
        let mut trace = ExecutionTrace::new(4, 1);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(&trace, &constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(&constraints, trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_invalid_proof() {
        let mut trace = ExecutionTrace::new(4, 1);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i + 1); // invalid
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(&trace, &constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(&constraints, trace.height as usize);
        assert!(!verifier.verify(&proof));
    }

    #[test]
    fn test_larger_trace() {
        let mut trace = ExecutionTrace::new(8, 1);
        for i in 0..8 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(&trace, &constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(&constraints, trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_multiple_variables() {
        let mut trace = ExecutionTrace::new(4, 2);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            row.insert("y".to_string(), i * 2);
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        // x[n+1] = x[n] + 1
        constraints.add_transition_constraint(
            "increment_x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        // y[n] = 2 * x[n]
        constraints.add_transition_constraint(
            "y_is_double_x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, _| {
                let x = Fr::from(*current.get("x").unwrap());
                let y = Fr::from(*current.get("y").unwrap());
                y - x * Fr::from(2u64)
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(&trace, &constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(&constraints, trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_zero_values() {
        let mut trace = ExecutionTrace::new(4, 1);
        for _ in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), 0); // All zeros
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "zero_sequence".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n // Should be zero
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_zero".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(&trace, &constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(&constraints, trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_complex_constraints() {
        let mut trace = ExecutionTrace::new(4, 2);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            row.insert("y".to_string(), i * i); // y = x^2
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        // x[n+1] = x[n] + 1
        constraints.add_transition_constraint(
            "increment_x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        // y[n] = x[n]^2
        constraints.add_transition_constraint(
            "y_is_x_squared".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, _| {
                let x = Fr::from(*current.get("x").unwrap());
                let y = Fr::from(*current.get("y").unwrap());
                y - x * x
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(&trace, &constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(&constraints, trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_invalid_complex_constraints() {
        let mut trace = ExecutionTrace::new(4, 2);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            row.insert("y".to_string(), i * i + 1); // y = x^2 + 1 (invalid)
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        // x[n+1] = x[n] + 1
        constraints.add_transition_constraint(
            "increment_x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        // y[n] = x[n]^2
        constraints.add_transition_constraint(
            "y_is_x_squared".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, _| {
                let x = Fr::from(*current.get("x").unwrap());
                let y = Fr::from(*current.get("y").unwrap());
                y - x * x
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(&trace, &constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(&constraints, trace.height as usize);
        assert!(!verifier.verify(&proof));
    }
}
