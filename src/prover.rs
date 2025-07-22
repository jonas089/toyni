use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::merkle::MerkleTree;
use crate::{digest_sha2, program::trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use num_bigint::BigUint;
use num_traits::ToPrimitive;

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
}

impl StarkProver {
    /// Creates a new STARK prover for the given trace and constraints.
    ///
    /// # Arguments
    ///
    /// * `trace` - The execution trace to prove
    /// * `constraints` - The constraint system defining program rules
    pub fn new(trace: ExecutionTrace) -> Self {
        Self { trace }
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
        let trace_len = self.trace.trace.len() as usize;
        let domain = GeneralEvaluationDomain::<Fr>::new(trace_len).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 8).unwrap();

        for column_idx in 0..self.trace.trace.len() {
            //let poly = self.trace.interpolate_column(domain, column_idx);
        }
        // todo: interpolate the trace columns one by one
        // Ci(x) = T(gx) - T(x) over the original domain for each T (note this is just an example constraint, the system must be adjusted accordingly to handle different kinds of constraints)
        // C(x) = C1(x) + C2(x) + ... + Cn(x)
        // Note that Ci(x) is always the same for the given program, and it's a symbolic polynomial meaning it does not hold any logic from the verifier point of view.
        // todo: implement proper air table

        // Create vanishing polynomial
        /*let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());

        // Divide to get quotient polynomial
        let (quotient_poly, _) = c_poly.divide(&z_poly).unwrap();

        // Evaluate quotient polynomial over extended domain
        let mut q_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| quotient_poly.evaluate(x))
            .collect();

        // todo: take each individual constraint polynomial and evaluate them over the extended domain,
        // committing the results to a merkle tree.
        // Share the functions that were used to interpolate the constraint polynomials with the verifier.
        // Sample random points with fiat shamir and check that C(z) = C'(z) in the verifier code.

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

        // Generate random challenges for verification
        let proof_transcript = build_proof_transcript(
            &q_evals,
            &fri_layers,
            &fri_challenges,
            &c_poly,
            &folding_commitment_trees,
        );

        let verifier_random_challenges =
            generate_spot_check_challenges(&proof_transcript, &extended_domain, total_queries);

        StarkProof {
            quotient_eval_domain: fri_layers[0].clone(),
            fri_layers,
            fri_challenges,
            combined_constraint: c_poly,
            quotient_poly,
            folding_commitment_trees,
            verifier_random_challenges,
        }*/
        StarkProof {
            quotient_eval_domain: vec![],
            fri_layers: vec![],
            fri_challenges: vec![],
            combined_constraint: ToyniPolynomial::new(vec![Fr::from(0)]),
            quotient_poly: ToyniPolynomial::new(vec![Fr::from(0)]),
            folding_commitment_trees: vec![],
            verifier_random_challenges: vec![],
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
