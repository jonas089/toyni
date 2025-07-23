use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::merkle::MerkleTree;
use crate::{digest_sha2, program::trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, BigInteger, Field, PrimeField, UniformRand};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain};
use rand::thread_rng;

fn random_poly(degree: usize) -> ToyniPolynomial {
    let mut rng = thread_rng();
    let coeffs: Vec<Fr> = (0..=degree).map(|_| Fr::rand(&mut rng)).collect();
    ToyniPolynomial::new(coeffs)
}

#[derive(Debug)]
pub struct StarkProof {
    pub quotient_eval_domain: Vec<Fr>,
    pub fri_layers: Vec<Vec<Fr>>,
    pub fri_challenges: Vec<Fr>,
    pub combined_constraint: ToyniPolynomial,
    pub quotient_poly: ToyniPolynomial,
    pub folding_commitment_trees: Vec<MerkleTree>,
    pub verifier_random_challenges: Vec<Fr>,
}

pub struct StarkProver {
    trace: ExecutionTrace,
}

impl StarkProver {
    pub fn new(trace: ExecutionTrace) -> Self {
        Self { trace }
    }

    pub fn generate_proof(&self) -> StarkProof {
        let trace_len = self.trace.trace.len() as usize;

        let domain = GeneralEvaluationDomain::<Fr>::new(trace_len).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 8).unwrap();
        let domain_slice: Vec<Fr> = domain.elements().collect();

        let mut trace_polys = Vec::new();

        for column_idx in 0..self.trace.trace[0].len() {
            let poly = self.trace.interpolate_column(&domain_slice, column_idx);
            trace_polys.push(poly);
        }

        fn fibonacci_constraint(ti2: Fr, ti1: Fr, ti0: Fr) -> Fr {
            ti2 - (ti1 + ti0)
        }

        let mut domain_evaluations = vec![Fr::ZERO; domain.size() - 2];

        for i in 0..(domain.size() - 2) {
            let ti2 = trace_polys[0](domain.element(i + 2));
            let ti1 = trace_polys[0](domain.element(i + 1));
            let ti0 = trace_polys[0](domain.element(i));
            domain_evaluations[i] = fibonacci_constraint(ti2, ti1, ti0);
        }

        let ci_poly = ToyniPolynomial::from_dense_poly(DensePolynomial::from_coefficients_slice(
            &domain.ifft(&domain_evaluations),
        ));

        let transcript_seed = digest_sha2(
            &ci_poly
                .coefficients()
                .iter()
                .flat_map(|c| c.into_bigint().to_bytes_be())
                .collect::<Vec<_>>(),
        );
        let mut alpha_bytes = [0u8; 32];
        alpha_bytes.copy_from_slice(&transcript_seed[..32]);
        //let alpha = Fr::from_le_bytes_mod_order(&alpha_bytes);

        let c_poly = ci_poly; //.scale(alpha);

        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());

        let r_poly = random_poly(2);
        let c_z_poly = c_poly.add(&r_poly.mul(&z_poly));

        let (quotient_poly, _) = c_z_poly.divide(&z_poly).unwrap();

        let mut q_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| quotient_poly.evaluate(x))
            .collect();

        let mut fri_layers = vec![q_evals.clone()];
        let mut fri_challenges = Vec::new();
        let mut folding_commitment_trees: Vec<MerkleTree> = Vec::new();
        let mut fri_transcript = Vec::new();

        while q_evals.len() > 8 {
            for eval in &q_evals {
                fri_transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
            }

            let fri_hash = digest_sha2(&fri_transcript);
            let mut beta_bytes = [0u8; 32];
            beta_bytes.copy_from_slice(&fri_hash[..32]);
            let beta = Fr::from_le_bytes_mod_order(&beta_bytes);
            fri_challenges.push(beta);

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

        StarkProof {
            quotient_eval_domain: vec![],
            fri_layers,
            fri_challenges,
            combined_constraint: c_z_poly,
            quotient_poly,
            folding_commitment_trees,
            verifier_random_challenges: vec![],
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{program::trace::ExecutionTrace, prover::StarkProver, verifier::StarkVerifier};
    use ark_bls12_381::Fr;

    #[test]
    fn test_constraint_poly() {
        let mut execution_trace = ExecutionTrace::new();
        execution_trace.insert_column(vec![
            Fr::from(1),
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(5),
            Fr::from(8),
            Fr::from(13),
            Fr::from(21),
        ]);
        let stark = StarkProver::new(execution_trace.clone());
        let proof = stark.generate_proof();
        let verifier = StarkVerifier::new(execution_trace.trace.len());
        assert!(verifier.verify(&proof));
    }
}
