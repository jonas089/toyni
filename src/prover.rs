use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::merkle::MerkleTree;
use crate::{digest_sha2, program::trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, BigInteger, PrimeField, UniformRand};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial};
use rand::thread_rng;

#[allow(dead_code)]
fn random_poly(degree: usize) -> ToyniPolynomial {
    let mut rng = thread_rng();
    let coeffs: Vec<Fr> = (0..=degree).map(|_| Fr::rand(&mut rng)).collect();
    ToyniPolynomial::new(coeffs)
}

pub const CONSTRAINT_SPOT_CHECKS: usize = 8;

#[derive(Debug)]
pub struct StarkProof {
    pub fri_layers: Vec<Vec<Fr>>,
    pub fri_challenges: Vec<Fr>,
    pub quotient_poly: ToyniPolynomial,
    pub folding_commitment_trees: Vec<MerkleTree>,
    pub trace_spot_checks: [[Fr; 3]; CONSTRAINT_SPOT_CHECKS],
    pub constraint_spot_checks: [Fr; CONSTRAINT_SPOT_CHECKS],
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

        let trace_poly = self
            .trace
            .interpolate_column(&domain.elements().collect::<Vec<Fr>>(), 0);
        let extended_points = extended_domain.elements().collect::<Vec<Fr>>();

        /*let trace_lde = extended_points
        .iter()
        .map(|x| trace_poly.evaluate(*x))
        .collect::<Vec<Fr>>();*/

        let z_poly = domain.vanishing_polynomial();

        let g = extended_domain.group_gen();

        let quotient_points: Vec<Fr> = extended_points
            .iter()
            .map(|&w| {
                let z = z_poly.evaluate(&w);
                if z == Fr::ZERO {
                    Fr::ZERO // or panic/assert if this should never happen
                } else {
                    let t0 = trace_poly.evaluate(w);
                    let t1 = trace_poly.evaluate(g * w);
                    let t2 = trace_poly.evaluate(g * g * w);
                    fibonacci_constraint(t2, t1, t0) / z
                }
            })
            .collect();

        let quotient_poly =
            DensePolynomial::from_coefficients_slice(&extended_domain.ifft(&quotient_points));

        let mut q_evals = quotient_points.clone();

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

        // spot check the first N points
        // todo: build a merkle tree from the evaluations and use fiat shamir
        // to reveal part of the trace without a clear context / position
        let mut trace_spot_checks = [[Fr::ZERO; 3]; CONSTRAINT_SPOT_CHECKS];
        let mut constraint_spot_checks = [Fr::ZERO; CONSTRAINT_SPOT_CHECKS];
        for i in 0..CONSTRAINT_SPOT_CHECKS {
            let t0 = trace_poly.evaluate(extended_domain.element(i));
            let t1 = trace_poly.evaluate(extended_domain.element(i + 1));
            let t2 = trace_poly.evaluate(extended_domain.element(i + 2));
            trace_spot_checks[i] = [t0, t1, t2];

            constraint_spot_checks[i] = fibonacci_constraint(t2, t1, t0)
        }

        StarkProof {
            fri_layers,
            fri_challenges,
            quotient_poly: ToyniPolynomial::from_dense_poly(quotient_poly),
            folding_commitment_trees,
            trace_spot_checks,
            constraint_spot_checks,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{program::trace::ExecutionTrace, prover::StarkProver, verifier::StarkVerifier};
    use ark_bls12_381::Fr;

    #[test]
    fn test_fibonacci() {
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
            Fr::from(34),
            Fr::from(55),
            Fr::from(89),
            Fr::from(144),
            Fr::from(233),
            Fr::from(377),
            Fr::from(610),
            Fr::from(987),
        ]);
        let stark = StarkProver::new(execution_trace.clone());
        let proof = stark.generate_proof();
        let verifier = StarkVerifier::new(execution_trace.trace.len());
        assert!(verifier.verify(&proof));
    }
}
