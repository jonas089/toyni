//! STARK (Scalable Transparent Argument of Knowledge) implementation.
//!
//! Converts program execution into trace polynomial, constructs composition polynomial,
//! and uses FRI protocol to prove low-degree.

use crate::digest_sha2;
use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::merkle::MerkleTree;
use crate::vm::{constraints::ConstraintSystem, trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_poly::DenseUVPolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::thread_rng;

const VERIFIER_QUERIES: usize = 80;

/// STARK proof containing components for verification.
pub struct StarkProof {
    /// Quotient polynomial evaluations over domain
    pub quotient_eval_domain: Vec<Fr>,
    /// FRI protocol layers with folded evaluations
    pub fri_layers: Vec<Vec<Fr>>,
    /// Random challenges for FRI folding
    pub fri_challenges: Vec<Fr>,
    /// Combined constraint polynomial
    pub combined_constraint: ToyniPolynomial,
    /// Quotient polynomial from division
    pub quotient_poly: ToyniPolynomial,
    /// Merkle tree for quotient polynomial evaluations
    pub folding_commitment_trees: Vec<MerkleTree>,
    /// Fiat-Shamir random challenges
    pub verifier_random_challenges: Vec<Fr>,
}

/// STARK prover component.
pub struct StarkProver<'a> {
    /// Execution trace to prove
    trace: &'a ExecutionTrace,
    /// Constraint system defining rules
    constraints: &'a ConstraintSystem,
}

impl<'a> StarkProver<'a> {
    /// Creates new STARK prover.
    pub fn new(trace: &'a ExecutionTrace, constraints: &'a ConstraintSystem) -> Self {
        Self { trace, constraints }
    }

    /// Generates STARK proof for execution trace.
    pub fn generate_proof(&self) -> StarkProof {
        let trace_len = self.trace.height as usize;
        let domain = GeneralEvaluationDomain::<Fr>::new(trace_len).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 2).unwrap();

        println!("\n=== Prover Debug ===");
        println!("Trace length: {}", trace_len);
        println!("Extended domain size: {}", extended_domain.size());
        let constraint_polys = self.constraints.interpolate_all_constraints(self.trace);
        println!("\nConstraint polynomials:");
        for (i, poly) in constraint_polys.iter().enumerate() {
            println!("Constraint {}: {:?}", i, poly.coefficients);
        }

        let combined_constraint = constraint_polys
            .iter()
            .fold(ToyniPolynomial::zero(), |acc, p| acc.add(p));
        println!(
            "\nCombined constraint: {:?}",
            combined_constraint.coefficients
        );

        let c_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| combined_constraint.evaluate(x))
            .collect();
        println!("\nConstraint evaluations at domain points:");
        for (i, eval) in c_evals.iter().enumerate() {
            println!("C[{}] = {:?}", i, eval);
        }

        let c_poly = DensePolynomial::from_coefficients_slice(&extended_domain.ifft(&c_evals));
        let c_poly = ToyniPolynomial::from_dense_poly(c_poly);
        println!(
            "\nInterpolated constraint polynomial: {:?}",
            c_poly.coefficients
        );

        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());
        // todo: add a random polynomial to the vanishing polynomial
        println!("\nVanishing polynomial: {:?}", z_poly.coefficients);

        let (quotient_poly, rem) = c_poly.divide(&z_poly).unwrap();
        println!("\nQuotient polynomial: {:?}", quotient_poly.coefficients);
        println!("Remainder polynomial: {:?}", rem.coefficients);

        let mut q_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| quotient_poly.evaluate(x))
            .collect();
        println!("\nQuotient evaluations at domain points:");
        for (i, eval) in q_evals.iter().enumerate() {
            println!("Q[{}] = {:?}", i, eval);
        }

        let mut fri_layers = vec![q_evals.clone()];
        let mut fri_challenges = Vec::new();
        // large enough degree to ensure enough points are checked
        let mut folding_commitment_trees: Vec<MerkleTree> = Vec::new();
        while q_evals.len() > 4 {
            let beta = Fr::rand(&mut thread_rng());
            fri_challenges.push(beta);
            // todo: commit the q_evals to the merkle tree
            q_evals = fri_fold(&q_evals, beta);
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

        // todo: actually hash the real proof data
        // this requires a serializable Fr type
        // currently this is a mock implementation of non-interactive fiat shamir
        // it must be improved to be secure and circuit-specific
        let proof_hash = digest_sha2(&[0; 32]);
        let proof_hash_u32: Vec<u32> = proof_hash
            .chunks(4)
            .map(|chunk| {
                let mut buf = [0u8; 4];
                buf.copy_from_slice(chunk);
                u32::from_be_bytes(buf)
            })
            .collect(); // ✅ now collects into Vec<u32>
        let proof_hash_int = BigUint::from_slice(&proof_hash_u32);
        let mut verifier_random_challenges = Vec::new();
        for _ in 0..VERIFIER_QUERIES {
            let verifier_index = proof_hash_int.clone() + BigUint::from(1 as u64);
            let verifier_index_field = verifier_index % BigUint::from(extended_domain.size());
            let random_interactive_challenge =
                extended_domain.element(verifier_index_field.to_usize().unwrap());
            verifier_random_challenges.push(random_interactive_challenge);
        }

        StarkProof {
            // todo: commit the quotient_eval_domain to the merkle tree
            // don't reveal all evaluations, use Fiat-Shamir transform
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

/// STARK verifier component.
pub struct StarkVerifier<'a> {
    /// Constraint system defining rules
    #[allow(unused)]
    constraints: &'a ConstraintSystem,
    /// Length of execution trace
    trace_len: usize,
}

impl<'a> StarkVerifier<'a> {
    /// Creates new STARK verifier.
    pub fn new(constraints: &'a ConstraintSystem, trace_len: usize) -> Self {
        Self {
            constraints,
            trace_len,
        }
    }

    /// Verifies STARK proof.
    pub fn verify(&self, proof: &StarkProof) -> bool {
        let domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len * 2).unwrap();
        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());

        println!("\n=== Verifier Debug ===");
        println!("Trace length: {}", self.trace_len);
        println!("Extended domain size: {}", extended_domain.size());

        for i in proof.verifier_random_challenges.iter() {
            // todo: use Fiat-Shamir transform to get the random interactive challenge
            let random_interactive_challenge =
                extended_domain.element(rand::random::<usize>() % extended_domain.size());
            let q_eval = proof.quotient_poly.evaluate(random_interactive_challenge);
            let z_eval = z_poly.evaluate(random_interactive_challenge);
            let c_eval = proof
                .combined_constraint
                .evaluate(random_interactive_challenge);

            // f_next(x) = (f(x) + f(-x) + β * (f(x) - f(-x))) / 2
            // f₁[j] == (f₀[j] + f₀[j + n/2] + βᵢ * (f₀[j] - f₀[j + n/2])) / 2
            // todo: check if the merkle commitments are consistent with the evaluations
            // perform the consistency check at the random interactive challenge for all layers in the tree

            // in non-interactive mode we want to use the same indices for the spot checks in the
            // merkle commitment check. To achieve this we must hash the proof data
            // e.g. transcript = hash(public_inputs || quotient_merkle_root || fri_merkle_roots || ...)
            // and use that hash to derive the random queries deterministically.

            todo!("Check that the merkle tree commitments are consistent with the evaluations");

            println!("\nSpot check {}:", i);
            println!("x₀ = {:?}", random_interactive_challenge);
            println!("Q(x₀) = {:?}", q_eval);
            println!("Z(x₀) = {:?}", z_eval);
            println!("C(x₀) = {:?}", c_eval);
            println!("Q(x₀) * Z(x₀) = {:?}", q_eval * z_eval);

            if q_eval * z_eval != c_eval {
                println!("❌ Spot check failed: Q(x₀)*Z(x₀) ≠ C(x₀)");
                return false;
            }
        }

        // FRI folding check
        // f_{i+1}(x^2) = (f_i(x) + f_i(-x) + βᵢ * (f_i(x) - f_i(-x))) / 2
        let mut current = &proof.quotient_eval_domain;
        for (i, beta) in proof.fri_challenges.iter().enumerate() {
            let folded = fri_fold(current, *beta);
            if proof.fri_layers.get(i + 1).expect("FRI layer should exist") != &folded {
                println!("❌ FRI folding failed at layer {}", i);
                return false;
            }
            current = &proof.fri_layers[i + 1];
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
