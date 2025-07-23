use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, BigInteger, Field, PrimeField};
use ark_poly::{
    DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial,
};

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

        // Hardcoded alpha for now
        //let alpha = Fr::from(7u64);

        // Recompute Fibonacci constraint ci(x) over the domain
        let elements: Vec<Fr> = domain.elements().collect();
        let mut ci_evals = vec![Fr::ZERO; domain.size()];

        for i in 0..(domain.size() - 2) {
            let x0 = elements[i];
            let x1 = elements[i + 1];
            let x2 = elements[i + 2];

            let a0 = proof.combined_constraint.evaluate(x0);
            let a1 = proof.combined_constraint.evaluate(x1);
            let a2 = proof.combined_constraint.evaluate(x2);

            // ci(x) = a2 - (a1 + a0)
            ci_evals[i] = a2 - (a1 + a0);
        }

        let ci_poly = Polynomial::from_dense_poly(DensePolynomial::from_coefficients_slice(
            &domain.ifft(&ci_evals),
        ));

        let expected_c_poly = ci_poly; //.scale(alpha);

        for i in 0..5 {
            let x = elements[i];
            let expected = expected_c_poly.evaluate(x);
            let actual = proof.combined_constraint.evaluate(x);

            println!("expected: {:?}, actual: {:?}", expected, actual);

            if expected != actual {
                println!(
                    "❌ combined_constraint(x) != α * ci(x) at i={}: {} != {}",
                    i, actual, expected
                );
                return false;
            }
        }

        // Fiat-Shamir random point checks: Q(x)·Z(x) == C(x)
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
            let x = Fr::from_le_bytes_mod_order(&challenge_bytes);

            let q_eval = proof.quotient_poly.evaluate(x);
            let z_eval = z_poly.evaluate(x);
            let c_eval = proof.combined_constraint.evaluate(x);

            if q_eval * z_eval != c_eval {
                println!("❌ Spot check failed: Q(x₀)*Z(x₀) ≠ C(x₀)");
                return false;
            }
        }

        if !self.verify_fri_layers(
            &proof.fri_layers,
            &proof.fri_challenges,
            &proof.folding_commitment_trees,
        ) {
            println!("❌ FRI folding consistency check failed");
            return false;
        }

        true
    }

    fn verify_fri_layers(
        &self,
        fri_layers: &[Vec<Fr>],
        fri_challenges: &[Fr],
        folding_commitment_trees: &[crate::merkle::MerkleTree],
    ) -> bool {
        let mut current_layer = &fri_layers[0];
        for (i, ((next_layer, beta), merkle_tree)) in fri_layers[1..]
            .iter()
            .zip(fri_challenges)
            .zip(folding_commitment_trees)
            .enumerate()
        {
            let half_n = current_layer.len() / 2;

            // maximum domain size is 64
            assert!(current_layer.len() < 65);

            if next_layer.len() != half_n {
                println!(
                    "❌ FRI layer {} has incorrect size: expected {}, got {}",
                    i,
                    half_n,
                    next_layer.len()
                );
                return false;
            }

            for j in 0..half_n {
                let fx = current_layer[j];
                let f_neg_x = current_layer[j + half_n];

                // Fold: (f(x) + f(-x) + β(f(x) - f(-x))) / 2
                let numerator = fx + f_neg_x + *beta * (fx - f_neg_x);
                let expected_next = numerator * Fr::from(2u64).inverse().unwrap();
                let actual_next = next_layer[j];

                if expected_next != actual_next {
                    println!("❌ FRI folding failed at layer {}, index {}", i, j);
                    println!("Expected: {:?}", expected_next);
                    println!("Actual:   {:?}", actual_next);
                    return false;
                }

                // Merkle commitment verification
                let value_bytes = actual_next.into_bigint().to_bytes_be();
                let proof = merkle_tree.get_proof(j).expect("Merkle proof should exist");
                let root = merkle_tree.root().expect("Merkle root should exist");

                if !crate::merkle::verify_merkle_proof(value_bytes, &proof, &root) {
                    println!(
                        "❌ Merkle proof verification failed at layer {}, index {}",
                        i, j
                    );
                    return false;
                }
            }

            current_layer = next_layer;
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
