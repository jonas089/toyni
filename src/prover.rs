use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::merkle::MerkleTree;
use crate::{digest_sha2, program::trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, BigInteger, PrimeField, UniformRand};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use rand::thread_rng;

#[allow(dead_code)]
fn random_poly(degree: usize) -> ToyniPolynomial {
    let mut rng = thread_rng();
    let coeffs: Vec<Fr> = (0..=degree).map(|_| Fr::rand(&mut rng)).collect();
    ToyniPolynomial::new(coeffs)
}

pub const CONSTRAINT_SPOT_CHECKS: usize = 50;

#[derive(Debug)]
pub struct StarkProof {
    pub fri_layers: Vec<Vec<Fr>>,
    pub fri_challenges: Vec<Fr>,
    pub deep_poly: ToyniPolynomial,
    pub folding_commitment_trees: Vec<MerkleTree>,
    pub trace_spot_checks: [[Fr; 3]; CONSTRAINT_SPOT_CHECKS],
    pub constraint_spot_checks: [Fr; CONSTRAINT_SPOT_CHECKS],
    pub r_poly: ToyniPolynomial,
    pub z: Fr,
    pub c_at_z: Fr,
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
        let shift = Fr::from(7);
        let extended_points = extended_domain.clone().elements().collect::<Vec<Fr>>();

        let r_poly = random_poly(2);

        fn fibonacci_constraint(t2: Fr, t1: Fr, t0: Fr) -> Fr {
            t2 - (t1 + t0)
        }

        let trace_poly = self
            .trace
            .interpolate_column(&domain.elements().collect::<Vec<Fr>>(), 0);

        let z_poly = domain.vanishing_polynomial();
        let z_toyni = ToyniPolynomial::from_dense_poly(z_poly.clone().into());
        let g = domain.group_gen();

        let fibonacci_column = self.trace.get_column(0);
        println!("fib column: {:?}", &fibonacci_column);

        // todo: don't interpolate the C poly, evaluate at points instead.
        // interpolation of valid constraints always yields 0
        let constraint_evals: Vec<Fr> = domain
            .elements()
            .take(domain.size() - 2)
            .map(|x| {
                let t0 = trace_poly.evaluate(x);
                let t1 = trace_poly.evaluate(g * x);
                let t2 = trace_poly.evaluate(g * g * x);
                fibonacci_constraint(t2, t1, t0)
            })
            .collect();

        let mut padded_constraints = constraint_evals.clone();
        while padded_constraints.len() < domain.size() {
            padded_constraints.push(Fr::ZERO);
        }

        let constraint_poly =
            Evaluations::from_vec_and_domain(padded_constraints, domain).interpolate_by_ref();

        let constraint_poly = ToyniPolynomial::from_dense_poly(constraint_poly);

        let z = Fr::rand(&mut thread_rng());
        let c_at_z = constraint_poly.evaluate(z);

        let deep_poly = constraint_poly
            .sub(&ToyniPolynomial::new(vec![c_at_z]))
            .divide_by_linear(z)
            .0
            .add(&r_poly.mul(&z_toyni));

        let deep_points: Vec<Fr> = extended_points
            .iter()
            .map(|&w| deep_poly.evaluate(w * shift))
            .collect();

        let mut d_evals = deep_points.clone();
        let mut fri_layers = vec![d_evals.clone()];
        let mut fri_challenges = Vec::new();
        let mut folding_commitment_trees: Vec<MerkleTree> = Vec::new();
        let mut fri_transcript = Vec::new();

        while d_evals.len() > 8 {
            for eval in &d_evals {
                fri_transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
            }

            let fri_hash = digest_sha2(&fri_transcript);
            let mut beta_bytes = [0u8; 32];
            beta_bytes.copy_from_slice(&fri_hash[..32]);
            let beta = Fr::from_le_bytes_mod_order(&beta_bytes);
            fri_challenges.push(beta);

            d_evals = fri_fold(&d_evals, beta);

            let folding_step_merkle_tree = MerkleTree::new(
                d_evals
                    .clone()
                    .iter()
                    .map(|x| x.into_bigint().to_bytes_be())
                    .collect(),
            );
            folding_commitment_trees.push(folding_step_merkle_tree);
            fri_layers.push(d_evals.clone());
        }

        let mut trace_spot_checks = [[Fr::ZERO; 3]; CONSTRAINT_SPOT_CHECKS];
        let mut constraint_spot_checks = [Fr::ZERO; CONSTRAINT_SPOT_CHECKS];
        for i in 0..CONSTRAINT_SPOT_CHECKS {
            let x = extended_domain.element(domain.size() + i) * shift;
            let t0 = trace_poly.evaluate(x);
            let t1 = trace_poly.evaluate(g * x);
            let t2 = trace_poly.evaluate(g * g * x);
            trace_spot_checks[i] = [t0, t1, t2];
            constraint_spot_checks[i] = fibonacci_constraint(t2, t1, t0);
        }

        println!("Trace Degree: {:?}", &trace_poly.degree());
        println!("Deep Degree: {:?}", &deep_poly.degree());
        println!("Constraint Degree: {:?}", &constraint_poly.degree());

        StarkProof {
            fri_layers,
            fri_challenges,
            deep_poly,
            folding_commitment_trees,
            trace_spot_checks,
            constraint_spot_checks,
            r_poly,
            z,
            c_at_z,
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
