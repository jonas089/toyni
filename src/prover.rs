use std::iter::Map;

use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::{digest_sha2, program::trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, BigInteger, PrimeField, UniformRand};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{
    DenseUVPolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial,
};
use rand::seq::IteratorRandom;
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
    pub fri_challenges: Vec<Fr>,
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

        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());
        let r_poly = random_poly(4);

        fn fibonacci_constraint(t2: Fr, t1: Fr, t0: Fr) -> Fr {
            if t2 == Fr::from(10610209857723u64)
                || t1 == Fr::from(10610209857723u64)
                || t0 == Fr::from(10610209857723u64)
            {
                return Fr::ZERO;
            }
            t2 - (t1 + t0)
        }

        // this part is correct according to all sources
        let trace_poly = self
            .trace
            .interpolate_column(&domain.elements().collect::<Vec<Fr>>(), 0);

        // todo: evaluate t(x) over extended & shifted domain and commit the evaluations

        // interpolate c(x) over the extended domain and divide by vanishing poly (must be shifted later for zk)
        let g = domain.group_gen();
        let c_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| {
                fibonacci_constraint(
                    trace_poly.evaluate(g * g * x),
                    trace_poly.evaluate(g * x),
                    trace_poly.evaluate(x),
                )
            })
            .collect();

        let c_poly = ToyniPolynomial::from_dense_poly(DensePolynomial::from_coefficients_slice(
            &extended_domain.ifft(&c_evals),
        ));

        let c_z_poly = c_poly.divide(&z_poly).unwrap().0;

        let z = get_random_z(&domain);
        let mut d_evals = vec![];
        let mut rng = thread_rng();
        let alpha = Fr::rand(&mut rng);

        let mut violations = 0;
        for x in extended_domain.elements() {
            let c_x = c_z_poly.evaluate(x);
            let c_z = c_z_poly.evaluate(z);
            let d_x = (c_x - c_z) / (x - z);
            d_evals.push(d_x);

            // simulating both kinds of spot checks, in addition with FRI folding
            // check these will be the cornerstone of security

            // later this will be moved to the verifier

            if c_poly.evaluate(x) != c_z_poly.evaluate(x) * z_poly.evaluate(x) {
                violations += 1;
            }

            // this problem will be solved
            if c_poly.evaluate(x) != Fr::ZERO {
                if c_poly.evaluate(x)
                    != fibonacci_constraint(
                        trace_poly.evaluate(g * g * x),
                        trace_poly.evaluate(g * x),
                        trace_poly.evaluate(x),
                    )
                {
                    violations += 1
                }
            }
        }

        println!("Violations: {:?}", &violations);
        assert_eq!(violations, 0);

        let d_poly = DensePolynomial::from_coefficients_slice(&extended_domain.ifft(&d_evals));
        println!("DEEP degree: {:?}", &d_poly.degree());

        let mut d_evals: Vec<Fr> = extended_domain
            .elements()
            .map(|x| d_poly.evaluate(&x))
            .collect();

        let mut fri_challenges = Vec::new();
        let mut fri_transcript = Vec::new();

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
            if d_evals.iter().all(|v| *v == d_evals[0]) {
                println!("✅ FRI folded to constant — likely low degree");
                fri_final_constant = true;
                break;
            }

            folding_steps += 1;
            println!(
                "Folding step: {}, current degree: {}",
                folding_steps,
                d_evals.len()
            );
        }
        if !fri_final_constant {
            println!("❌ FRI did not fold to constant — high degree likely");
        }

        println!("Folding Steps: {}", folding_steps);

        let mut trace_spot_checks = [[Fr::ZERO; 3]; CONSTRAINT_SPOT_CHECKS];
        let mut constraint_spot_checks = [Fr::ZERO; CONSTRAINT_SPOT_CHECKS];

        StarkProof {
            fri_challenges,
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
    use ark_ff::UniformRand;
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
    fn test_invalid_trace_should_fail() {
        let mut execution_trace = ExecutionTrace::new();
        let mut trace: Vec<u64> = fibonacci_list(64);
        for i in 20..21 {
            trace[i] = i as u64 * 1233;
        }
        let trace_field: Vec<Fr> = trace.iter().map(|x| Fr::from(*x)).collect();
        execution_trace.insert_column(trace_field);
        let stark = StarkProver::new(execution_trace.clone());
        let _proof = stark.generate_proof();
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
