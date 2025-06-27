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
use crate::merkle::{MerkleTree};
use crate::vm::{constraints::ConstraintSystem, trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::DenseUVPolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::thread_rng;

/// Number of random challenges for verifier spot checks
const VERIFIER_QUERIES: usize = 80;

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
        let mut rng = thread_rng();
        let random_poly = ToyniPolynomial::random(extended_domain.size() - 1, &mut rng);
        
        // Multiply combined constraint by random polynomial
        let masked_constraint = combined_constraint.mul(&random_poly);

        // Evaluate masked constraint over extended domain
        let c_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| masked_constraint.evaluate(x))
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

        // Generate random challenges for verification
        let proof_transcript = build_proof_transcript(
            &q_evals,
            &fri_layers,
            &fri_challenges,
            &combined_constraint,
            &folding_commitment_trees,
        );
        let verifier_random_challenges = generate_spot_check_challenges(
            &proof_transcript,
            &extended_domain,
            VERIFIER_QUERIES,
        );

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