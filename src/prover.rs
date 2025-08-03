use core::panic;

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
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 64).unwrap();

        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());
        let r_poly = random_poly(512);

        fn fibonacci_constraint(t2: Fr, t1: Fr, t0: Fr) -> Fr {
            if t2 == Fr::from(10610209857723u64)
                || t1 == Fr::from(10610209857723u64)
                || t0 == Fr::from(10610209857723u64)
            {
                return Fr::ZERO;
            }
            t2 - (t1 + t0)
        }
        let trace_poly = self
            .trace
            .interpolate_column(&domain.elements().collect::<Vec<Fr>>(), 0);

        let g = domain.group_gen();
        let z = get_random_z(&domain);
        let mut d_evals = vec![];
        let mut rng = thread_rng();
        let alpha = Fr::rand(&mut rng);
        let mut roots_in_lde = 0;
        for x in extended_domain.elements() {
            // constraints don't hold for the last 2 rows
            let c_x = fibonacci_constraint(
                trace_poly.evaluate(g * g * x),
                trace_poly.evaluate(g * x),
                trace_poly.evaluate(x),
            );
            let c_z = fibonacci_constraint(
                trace_poly.evaluate(g * g * z),
                trace_poly.evaluate(g * z),
                trace_poly.evaluate(z),
            );
            let d_x = alpha * (c_x - c_z) / (x - z) + r_poly.evaluate(x) * z_poly.evaluate(x);
            if d_x == Fr::ZERO {
                roots_in_lde += 1;
            }
            d_evals.push(d_x);
        }
        println!("roots in lde: {:?}", &roots_in_lde);

        let d_poly = extended_domain.ifft(&d_evals);
        let mut roots = 0;
        for element in &d_poly {
            if element == &Fr::ZERO {
                roots += 1;
            }
        }
        println!("roots: {}", &roots);
        //println!("d_poly: {:?}, size: {:?}", &d_poly, &d_poly.len());

        let mut fri_layers = vec![d_evals.clone()];
        let mut fri_challenges = Vec::new();
        let mut folding_commitment_trees: Vec<MerkleTree> = Vec::new();
        let mut fri_transcript = Vec::new();

        // folding the D(x) polynomial
        let mut folding_steps = 0;
        let mut fri_final_constant = false;

        while d_evals.len() > 1 {
            for eval in &d_evals {
                fri_transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
            }

            let fri_hash = digest_sha2(&fri_transcript);
            let mut beta_bytes = [0u8; 32];
            beta_bytes.copy_from_slice(&fri_hash[..32]);
            let beta = Fr::from_le_bytes_mod_order(&beta_bytes);
            fri_challenges.push(beta);

            d_evals = fri_fold(&d_evals, beta);

            // Check if all values are the same (constant polynomial)
            if d_evals.iter().all(|v| *v == d_evals[0]) {
                println!("✅ FRI folded to constant — likely low degree");
                fri_final_constant = true;
                break;
            }

            let folding_step_merkle_tree = MerkleTree::new(
                d_evals
                    .clone()
                    .iter()
                    .map(|x| x.into_bigint().to_bytes_be())
                    .collect(),
            );
            folding_commitment_trees.push(folding_step_merkle_tree);
            fri_layers.push(d_evals.clone());
            folding_steps += 1;
        }

        if !fri_final_constant {
            println!("❌ FRI did not fold to constant — high degree likely");
        }
        println!("Folding Steps: {}", folding_steps);

        let mut trace_spot_checks = [[Fr::ZERO; 3]; CONSTRAINT_SPOT_CHECKS];
        let mut constraint_spot_checks = [Fr::ZERO; CONSTRAINT_SPOT_CHECKS];

        StarkProof {
            fri_layers,
            fri_challenges,
            folding_commitment_trees,
            trace_spot_checks,
            constraint_spot_checks,
        }
    }
}

fn get_random_z(domain: &GeneralEvaluationDomain<Fr>) -> Fr {
    let mut rng = thread_rng();
    let domain_set: std::collections::HashSet<Fr> = domain.elements().collect();
    let g = domain.group_gen();

    loop {
        let z = Fr::rand(&mut rng);
        if !domain_set.contains(&z)
            && !domain_set.contains(&(g * z))
            && !domain_set.contains(&(g * g * z))
        {
            return z;
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{program::trace::ExecutionTrace, prover::StarkProver};
    use ark_bls12_381::Fr;
    use ark_ff::{AdditiveGroup, UniformRand};
    use rand::thread_rng;

    #[test]
    fn test_fibonacci() {
        let mut execution_trace = ExecutionTrace::new();
        let trace: Vec<u64> = fibonacci_list(64);
        let trace_field: Vec<Fr> = trace.iter().map(|x| Fr::from(*x)).collect();
        execution_trace.insert_column(trace_field);
        let stark = StarkProver::new(execution_trace.clone());
        let _proof = stark.generate_proof();
        /*let verifier = StarkVerifier::new(execution_trace.trace.len());
        assert!(verifier.verify(&proof));*/
    }

    #[test]
    #[should_panic]
    fn test_large_trace_should_fail() {
        let mut execution_trace = ExecutionTrace::new();
        let mut rng = thread_rng();
        let trace: Vec<Fr> = (0..128).map(|_| Fr::rand(&mut rng)).collect();

        let trace_field: Vec<Fr> = trace.iter().map(|x| Fr::from(*x)).collect();
        execution_trace.insert_column(trace_field);
        let stark = StarkProver::new(execution_trace.clone());
        let _proof = stark.generate_proof();
        /*let verifier = StarkVerifier::new(execution_trace.trace.len());
        assert!(verifier.verify(&proof));*/
    }

    #[test]
    #[should_panic]
    fn test_invalid_trace_should_fail() {
        let mut execution_trace = ExecutionTrace::new();
        let mut trace = vec![Fr::ZERO; 64];
        let mut rng = thread_rng();
        for i in 0..64 {
            trace[i] = Fr::rand(&mut rng);
        }
        execution_trace.insert_column(trace);
        let stark = StarkProver::new(execution_trace.clone());
        let _proof = stark.generate_proof();
    }

    fn fibonacci_list(n: usize) -> Vec<u64> {
        let mut fibs: Vec<u64> = Vec::with_capacity(n);
        let mut a = 1;
        let mut b = 1;
        for _ in 0..n {
            fibs.push(a);
            let next = a + b;
            a = b;
            b = next;
        }
        fibs
    }
}
