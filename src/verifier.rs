use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::{math::polynomial::Polynomial, merkle::verify_merkle_proof, prover::StarkProof, digest_sha2};

/// STARK verifier component that verifies proofs.
///
/// The verifier:
/// 1. Checks FRI folding consistency with Merkle proofs
/// 2. Verifies constraint satisfaction at random points
/// 3. Ensures all commitments are valid
pub struct StarkVerifier {
    /// Length of execution trace
    trace_len: usize,
}

impl StarkVerifier {
    /// Creates a new STARK verifier for the given trace length.
    ///
    /// # Arguments
    ///
    /// * `trace_len` - The length of the execution trace
    pub fn new(trace_len: usize) -> Self {
        Self {
            trace_len,
        }
    }

    /// Verifies a STARK proof.
    ///
    /// The verification process:
    /// 1. Checks FRI folding consistency with Merkle proofs
    /// 2. Verifies constraint satisfaction at random points
    /// 3. Ensures all commitments are valid
    /// 4. Verifies low degree of the quotient polynomial
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
        let _extended_domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len * 8).unwrap();
        let z_poly = Polynomial::from_dense_poly(domain.vanishing_polynomial().into());

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

        // Low degree check: verify the final FRI layer has the expected small degree
        let final_layer = proof.fri_layers.last().expect("FRI layers should not be empty");
        let expected_final_degree = 3; // The final layer should have degree ≤ 3
        let actual_final_degree = final_layer.len() - 1;
        
        if actual_final_degree > expected_final_degree {
            println!("❌ Low degree check failed: final layer has degree {} > {}", 
                     actual_final_degree, expected_final_degree);
            return false;
        }

        // Verify constraint satisfaction at random points
        let verifier_transcript = build_verifier_transcript(
            &proof.quotient_eval_domain,
            &proof.fri_layers,
            &proof.fri_challenges,
            &proof.combined_constraint,
            &proof.folding_commitment_trees,
        );
        let verifier_hash = digest_sha2(&verifier_transcript);
        
        for i in 0..proof.verifier_random_challenges.len() {
            // Generate deterministic challenge using Fiat-Shamir
            let mut challenge_bytes = [0u8; 32];
            let hash_offset = (i * 32) % verifier_hash.len();
            let bytes_to_copy = std::cmp::min(32, verifier_hash.len() - hash_offset);
            challenge_bytes[..bytes_to_copy].copy_from_slice(&verifier_hash[hash_offset..hash_offset + bytes_to_copy]);
            let random_interactive_challenge = Fr::from_le_bytes_mod_order(&challenge_bytes);
            
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

/// Builds a transcript for Fiat-Shamir challenge generation (verifier side)
fn build_verifier_transcript(
    quotient_eval_domain: &[Fr],
    fri_layers: &[Vec<Fr>],
    fri_challenges: &[Fr],
    combined_constraint: &Polynomial,
    folding_commitment_trees: &[crate::merkle::MerkleTree],
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