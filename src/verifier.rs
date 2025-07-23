use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::{digest_sha2, math::polynomial::Polynomial, prover::StarkProof};

pub struct StarkVerifier {
    trace_len: usize,
}

impl StarkVerifier {
    pub fn new(trace_len: usize) -> Self {
        Self { trace_len }
    }

    pub fn verify(&self, proof: &StarkProof) -> bool {
        let domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len).unwrap();
        let z_poly = Polynomial::from_dense_poly(domain.vanishing_polynomial().into());

        let verifier_transcript = build_verifier_transcript(
            &proof.quotient_eval_domain,
            &proof.fri_layers,
            &proof.fri_challenges,
            &proof.combined_constraint,
            &proof.folding_commitment_trees,
        );
        let verifier_hash = digest_sha2(&verifier_transcript);

        for i in 0..64 {
            let mut challenge_bytes = [0u8; 32];
            let hash_offset = (i * 32) % verifier_hash.len();
            let bytes_to_copy = std::cmp::min(32, verifier_hash.len() - hash_offset);
            challenge_bytes[..bytes_to_copy]
                .copy_from_slice(&verifier_hash[hash_offset..hash_offset + bytes_to_copy]);
            let random_challenge = Fr::from_le_bytes_mod_order(&challenge_bytes);

            let q_eval = proof.quotient_poly.evaluate(random_challenge);
            let z_eval = z_poly.evaluate(random_challenge);
            let c_eval = proof.combined_constraint.evaluate(random_challenge);

            if q_eval * z_eval != c_eval {
                println!("❌ Spot check failed: Q(x₀)*Z(x₀) ≠ C(x₀)");
                return false;
            }
        }

        true
    }
}

fn build_verifier_transcript(
    quotient_eval_domain: &[Fr],
    fri_layers: &[Vec<Fr>],
    fri_challenges: &[Fr],
    combined_constraint: &Polynomial,
    folding_commitment_trees: &[crate::merkle::MerkleTree],
) -> Vec<u8> {
    let mut transcript = Vec::new();

    for eval in quotient_eval_domain {
        transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
    }

    for layer in fri_layers {
        for eval in layer {
            transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
        }
    }

    for challenge in fri_challenges {
        transcript.extend_from_slice(&challenge.into_bigint().to_bytes_be());
    }

    for coeff in combined_constraint.coefficients() {
        transcript.extend_from_slice(&coeff.into_bigint().to_bytes_be());
    }

    for tree in folding_commitment_trees {
        if let Some(root) = tree.root() {
            transcript.extend_from_slice(&root);
        }
    }

    transcript
}
