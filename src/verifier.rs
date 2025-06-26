use ark_bls12_381::Fr;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::{math::polynomial::Polynomial, merkle::verify_merkle_proof, prover::StarkProof, vm::constraints::ConstraintSystem};

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
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(self.trace_len * 8).unwrap();
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