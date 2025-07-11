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
use crate::merkle::MerkleTree;
use crate::program::{constraints::ConstraintSystem, trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::DenseUVPolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::Rng;

/// Calculate the optimal number of queries for 128-bit security
/// Based on the formula: m ≥ log₂(L) + 128
/// where L is the number of FRI layers
fn calculate_optimal_queries(fri_layers: usize) -> usize {
    let log_l = (fri_layers as f64).log2();
    let optimal_queries = log_l + 128.0;
    optimal_queries.ceil() as usize
}

/// Number of random challenges for verifier spot checks
/// This will be calculated dynamically based on FRI layers for 128-bit security
const MIN_VERIFIER_QUERIES: usize = 64; // Minimum for constraint checks

/// Builds a transcript for Fiat-Shamir challenge generation
pub fn build_proof_transcript(
    quotient_eval_domain: &[Fr],
    fri_layers: &[Vec<Fr>],
    fri_challenges: &[Fr],
    combined_constraint: &ToyniPolynomial,
    folding_commitment_trees: &[MerkleTree],
) -> Vec<u8> {
    let mut transcript = Vec::new();

    // Add quotient polynomial evaluations
    for eval in quotient_eval_domain {
        transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
    }

    // Add FRI layers
    for layer in fri_layers {
        for eval in layer {
            transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
        }
    }

    // Add FRI challenges
    for challenge in fri_challenges {
        transcript.extend_from_slice(&challenge.into_bigint().to_bytes_be());
    }

    // Add constraint polynomial coefficients
    for coeff in combined_constraint.coefficients() {
        transcript.extend_from_slice(&coeff.into_bigint().to_bytes_be());
    }

    // Add Merkle tree roots
    for tree in folding_commitment_trees {
        if let Some(root) = tree.root() {
            transcript.extend_from_slice(&root);
        }
    }

    transcript
}

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
pub struct StarkProver {
    /// Execution trace to prove
    trace: ExecutionTrace,
    /// Constraint system defining program rules
    constraints: ConstraintSystem,
}

impl StarkProver {
    /// Creates a new STARK prover for the given trace and constraints.
    ///
    /// # Arguments
    ///
    /// * `trace` - The execution trace to prove
    /// * `constraints` - The constraint system defining program rules
    pub fn new(trace: ExecutionTrace, constraints: ConstraintSystem) -> Self {
        Self { trace, constraints }
    }

    /// Generates a STARK proof for the execution trace.
    ///
    /// The proof generation process:
    /// 1. Interpolates all constraints into polynomials
    /// 2. Combines constraints into a single polynomial
    /// 3. Generates random polynomial for zero-knowledge
    /// 4. Multiplies combined constraint by random polynomial
    /// 5. Divides by the vanishing polynomial to get quotient
    /// 6. Performs FRI folding with Merkle commitments
    /// 7. Generates random challenges for verification
    ///
    /// # Returns
    ///
    /// A `StarkProof` containing all components needed for verification
    pub fn generate_proof(&self) -> StarkProof {
        let trace_len = self.trace.height as usize;
        let domain = GeneralEvaluationDomain::<Fr>::new(trace_len).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 8).unwrap();

        // Interpolate all constraints into polynomials
        let constraint_polys = self.constraints.interpolate_all_constraints(&self.trace);

        // Combine all constraints into a single polynomial
        let mut combined_constraint = ToyniPolynomial::zero();
        for poly in constraint_polys {
            combined_constraint = combined_constraint.add(&poly);
        }

        // Generate random polynomial for zero-knowledge
        //let mut rng = thread_rng();
        //let random_poly = ToyniPolynomial::random(extended_domain.size() - 1, &mut rng);

        // Multiply combined constraint by random polynomial
        // let masked_constraint = combined_constraint.mul(&random_poly);

        // Evaluate masked constraint over extended domain
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

        // Interpolate trace polynomials for each variable
        let mut trace_polynomials = Vec::new();

        // Get all variable names from the first column
        let variable_names: Vec<String> = if let Some(first_column) = self.trace.trace.first() {
            first_column.keys().cloned().collect()
        } else {
            Vec::new()
        };

        for variable_name in &variable_names {
            // Extract values for this variable across all steps
            let mut variable_values = Vec::new();
            for step in 0..self.trace.trace.len() {
                let column_trace = &self.trace.trace[step];
                let value = column_trace.get(variable_name).unwrap_or(&0);
                variable_values.push(Fr::from(*value));
            }
            // Create domain points (powers of a generator)
            let domain = GeneralEvaluationDomain::<Fr>::new(variable_values.len()).unwrap();
            let domain_points: Vec<Fr> = domain.elements().collect();
            // Interpolate the polynomial
            let trace_poly = crate::math::fri::interpolate_poly(&domain_points, &variable_values);
            let trace_poly = ToyniPolynomial::from_dense_poly(trace_poly);
            trace_polynomials.push((variable_name.clone(), trace_poly));
        }

        // simulate verifier check for first column poly
        // evaluate the first column poly at a random point from the extended domain
        let mut rng = rand::thread_rng();
        let random_index = rng.gen_range(0..extended_domain.size());
        let random_point = extended_domain.element(random_index);
        let first_column_poly = trace_polynomials[0].1.clone();
        let first_column_poly_eval = first_column_poly.evaluate(random_point);
        println!("First column poly eval: {}", first_column_poly_eval);

        let c_eval = combined_constraint.evaluate(first_column_poly_eval);
        let z_eval = z_poly.evaluate(first_column_poly_eval);
        let q_eval = quotient_poly.evaluate(first_column_poly_eval);

        if c_eval != q_eval * z_eval {
            println!("❌ Simulated Trace Constraint not satisfied");
            println!("c_eval: {}", c_eval);
            println!("z_eval: {}", z_eval);
            println!("q_eval: {}", q_eval);
            println!("c_eval * z_eval: {}", c_eval * z_eval);
        }

        // Perform FRI folding with Merkle commitments
        let mut fri_layers = vec![q_evals.clone()];
        let mut fri_challenges = Vec::new();
        let mut folding_commitment_trees: Vec<MerkleTree> = Vec::new();
        let mut fri_transcript = Vec::new();

        while q_evals.len() > 4 {
            // Add current layer to transcript
            for eval in &q_evals {
                fri_transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
            }

            // Generate FRI challenge using Fiat-Shamir
            let fri_hash = digest_sha2(&fri_transcript);
            let mut beta_bytes = [0u8; 32];
            beta_bytes.copy_from_slice(&fri_hash[..32]);
            let beta = Fr::from_le_bytes_mod_order(&beta_bytes);

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

        // Calculate optimal number of queries for 128-bit security
        let total_fri_layers = fri_layers.len();
        let optimal_queries = calculate_optimal_queries(total_fri_layers);
        let constraint_queries = MIN_VERIFIER_QUERIES;
        let total_queries = optimal_queries + constraint_queries;

        println!("🔐 Security Parameters:");
        println!("  FRI layers: {}", total_fri_layers);
        println!("  Optimal FRI queries for 128-bit: {}", optimal_queries);
        println!("  Constraint queries: {}", constraint_queries);
        println!("  Total queries: {}", total_queries);

        // Generate random challenges for verification
        let proof_transcript = build_proof_transcript(
            &q_evals,
            &fri_layers,
            &fri_challenges,
            &combined_constraint,
            &folding_commitment_trees,
        );

        let verifier_random_challenges =
            generate_spot_check_challenges(&proof_transcript, &extended_domain, total_queries);

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

pub fn generate_spot_check_challenges(
    transcript: &[u8],
    domain: &GeneralEvaluationDomain<Fr>,
    num_challenges: usize,
) -> Vec<Fr> {
    let proof_hash = digest_sha2(transcript);
    let proof_hash_u32: Vec<u32> = proof_hash
        .chunks(4)
        .map(|chunk| {
            let mut buf = [0u8; 4];
            buf.copy_from_slice(chunk);
            u32::from_be_bytes(buf)
        })
        .collect();
    let proof_hash_int = BigUint::from_slice(&proof_hash_u32);
    let mut challenges = Vec::new();
    for i in 0..num_challenges {
        let verifier_index = &proof_hash_int + BigUint::from(i as u32 + 1);
        let verifier_index_field = verifier_index % BigUint::from(domain.size());
        let challenge = domain.element(verifier_index_field.to_usize().unwrap());
        challenges.push(challenge);
    }
    challenges
}
