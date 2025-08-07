use crate::prover::StarkProof;
use ark_bls12_381::Fr;
use ark_ff::{BigInteger, Field, PrimeField};

pub struct StarkVerifier {
    trace_len: usize,
}

const EXTENDED_DOMAIN_SIZE: usize = 128;

impl StarkVerifier {
    pub fn new(trace_len: usize) -> Self {
        Self { trace_len }
    }

    pub fn verify(&self, proof: &StarkProof) -> bool {
        // todo: DEEP verification using D(x) = (a * (C(x) - C(z)) / x - z) + r(x) * Zh(x) and openings of T(x), T(gx), T(ggx)

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

                let numerator = fx + f_neg_x + *beta * (fx - f_neg_x);
                let expected_next = numerator * Fr::from(2u64).inverse().unwrap();
                let actual_next = next_layer[j];

                if expected_next != actual_next {
                    println!("❌ FRI folding failed at layer {}, index {}", i, j);
                    println!("Expected: {:?}", expected_next);
                    println!("Actual:   {:?}", actual_next);
                    return false;
                }

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
