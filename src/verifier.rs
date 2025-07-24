use crate::{
    digest_sha2,
    math::polynomial::Polynomial,
    prover::{CI_SPOT_CHECKS, StarkProof},
};
use ark_bls12_381::Fr;
use ark_ff::{BigInteger, Field, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

pub struct StarkVerifier {
    trace_len: usize,
}

const FOLDING_SPOT_CHECKS: usize = 64;
const EXTENDED_DOMAIN_SIZE: usize = 128;

impl StarkVerifier {
    pub fn new(trace_len: usize) -> Self {
        Self { trace_len }
    }

    pub fn verify(&self, proof: &StarkProof) -> bool {
        let domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len * 8).unwrap();
        let z_poly = Polynomial::from_dense_poly(domain.vanishing_polynomial().into());
        let c_poly = proof.combined_constraint.sub(&proof.r_poly.mul(&z_poly));
        // should query random points over extended domain
        // and assert ci(x) == fib(ggx, gx, x)
        fn fibonacci_constraint(ti2: Fr, ti1: Fr, ti0: Fr) -> Fr {
            ti2 - (ti1 + ti0)
        }
        let ci_poly = proof.constraint_polys.first().unwrap();

        // check ci_poly correctness against fibonacci function
        // currently we only check the first 8 points, but we should use fiat shamir for this
        // sadly this method currently only works for the orginal domain,
        // meaning we are leaking trace values :( - must find a way to fix this!

        for i in 0..CI_SPOT_CHECKS {
            let trace_at_spot = proof.trace_spot_checks[i];
            let ti0 = trace_at_spot[0];
            let ti1 = trace_at_spot[1];
            let ti2 = trace_at_spot[2];

            let expected = fibonacci_constraint(ti2, ti1, ti0);
            let actual = ci_poly.evaluate(domain.element(i));

            assert_eq!(expected, actual);
        }

        let mut ci_transcript = Vec::new();
        for tree in &proof.folding_commitment_trees {
            if let Some(root) = tree.root() {
                ci_transcript.extend_from_slice(&root);
            }
        }
        for coeff in proof.combined_constraint.coefficients() {
            ci_transcript.extend_from_slice(&coeff.into_bigint().to_bytes_be());
        }

        let seed = digest_sha2(&ci_transcript);
        let mut seed_bytes = [0u8; 32];
        seed_bytes.copy_from_slice(&seed[..32]);

        // Fiat-Shamir random point checks: Q(x)·Z(x) == C(x)
        let verifier_transcript = build_verifier_transcript(
            &proof.fri_layers,
            &proof.fri_challenges,
            &proof.combined_constraint,
            &proof.folding_commitment_trees,
        );
        let verifier_hash = digest_sha2(&verifier_transcript);

        for i in 0..FOLDING_SPOT_CHECKS {
            let mut challenge_bytes = [0u8; 32];
            let hash_offset = (i * 32) % verifier_hash.len();
            let bytes_to_copy = std::cmp::min(32, verifier_hash.len() - hash_offset);
            challenge_bytes[..bytes_to_copy]
                .copy_from_slice(&verifier_hash[hash_offset..hash_offset + bytes_to_copy]);
            let x = Fr::from_le_bytes_mod_order(&challenge_bytes);

            let q_eval = proof.quotient_poly.evaluate(x);
            let z_eval = z_poly.evaluate(x);
            let c_eval = c_poly.evaluate(x);

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

            // maximum domain size is 128
            assert!(current_layer.len() < EXTENDED_DOMAIN_SIZE + 1);

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
    fri_layers: &[Vec<Fr>],
    fri_challenges: &[Fr],
    combined_constraint: &Polynomial,
    folding_commitment_trees: &[crate::merkle::MerkleTree],
) -> Vec<u8> {
    let mut transcript = Vec::new();

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
