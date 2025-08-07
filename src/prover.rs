use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::{digest_sha2, program::trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, BigInteger, Field, PrimeField, UniformRand};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use rand::thread_rng;

#[allow(dead_code)]
fn random_poly(degree: usize) -> ToyniPolynomial {
    let mut rng = thread_rng();
    let coeffs: Vec<Fr> = (0..=degree).map(|_| Fr::rand(&mut rng)).collect();
    ToyniPolynomial::new(coeffs)
}

#[derive(Debug)]
pub struct StarkProof {
    pub fri_challenges: Vec<Fr>,
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
        // the extended domain that overlaps the original roots of unity domain
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 8).unwrap();
        let shifted_domain = extended_domain.get_coset(Fr::from(7)).unwrap();

        // the vanishing polynomial for our trace over the original domain
        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());

        // the fibonacci constraint used for proving
        fn fibonacci_constraint(t2: Fr, t1: Fr, t0: Fr) -> Fr {
            // temporary solution - this is insecure and proper boundaries are better
            t2 - (t1 + t0)
        }

        fn boundary_constraint_1(x: Fr, g: Fr, n: usize) -> Fr {
            x - g.pow([(n - 1) as u64])
        }

        fn boundary_constraint_2(x: Fr, g: Fr, n: usize) -> Fr {
            x - g.pow([(n - 2) as u64])
        }

        // the trace interpolated as a polynomial over the original domain
        let trace_poly = self
            .trace
            .interpolate_column(&domain.elements().collect::<Vec<Fr>>(), 0);

        let g = domain.group_gen();

        let z = get_random_z(&extended_domain, &shifted_domain);

        let c_evals: Vec<Fr> = shifted_domain
            .elements()
            .map(|x| {
                fibonacci_constraint(
                    trace_poly.evaluate(g * g * x),
                    trace_poly.evaluate(g * x),
                    trace_poly.evaluate(x),
                ) * boundary_constraint_1(x, g, trace_len)
                    * boundary_constraint_2(x, g, trace_len)
            })
            .collect();

        let c_poly = ToyniPolynomial::new(shifted_domain.ifft(&c_evals));

        // evaluations of the quotient polynomial at challenge points
        let mut q_evals: Vec<Fr> = Vec::new();
        for x in shifted_domain.elements() {
            q_evals.push(c_poly.evaluate(x) / z_poly.evaluate(x));
        }

        // interpolation of the quotient for development purposes
        let q_poly = ToyniPolynomial::new(shifted_domain.ifft(&q_evals));
        let mut test_spot_check: Vec<Fr> = Vec::new();

        let q_z = q_poly.evaluate(z);
        let t_z = trace_poly.evaluate(z);
        let t_gz = trace_poly.evaluate(g * z);
        let t_ggz = trace_poly.evaluate(g * g * z);
        let mut d_evals = vec![];

        // evaluate DEEP polynomial
        for x in shifted_domain.elements() {
            let q_x = q_poly.evaluate(x);
            let t_x = trace_poly.evaluate(x);
            let t_gx = trace_poly.evaluate(g * x);
            let t_ggx = trace_poly.evaluate(g * g * x);
            let d_x = (q_x - q_z) / (x - z)
                + (t_ggx - t_ggz) / (x - z)
                + (t_gx - t_gz) / (x - z)
                + (t_x - t_z) / (x - z);
            test_spot_check.push(d_x.clone());
            d_evals.push(d_x);
        }

        // simulate consistency check at z
        assert_eq!(
            fibonacci_constraint(
                trace_poly.evaluate(g * g * z),
                trace_poly.evaluate(g * z),
                trace_poly.evaluate(z),
            ) * boundary_constraint_1(z, g, trace_len)
                * boundary_constraint_2(z, g, trace_len),
            q_poly.evaluate(z) * z_poly.evaluate(z)
        );

        // we fold the polynomial using our FRI evaluation domain
        // the spot checks will later ensure that the polynomial was folded correctly
        // and that it is of low degree
        let fri_challenges = Vec::new();
        let mut fri_transcript: Vec<u8> = Vec::new();
        let mut folding_steps: usize = 0;

        // note the degree of the deep polynomial for debugging
        let d_poly_degree: usize = ToyniPolynomial::new(shifted_domain.ifft(&d_evals)).degree();

        let mut xs: Vec<Fr> = shifted_domain.elements().collect();

        while d_evals.len() > 1 {
            // todo: merkle commit to all evaluations so verifier can check consistency
            // also verifier will compare d(x) to the spot check x mentioned above
            for eval in &d_evals {
                fri_transcript.extend_from_slice(&eval.into_bigint().to_bytes_be());
            }
            let fri_hash = digest_sha2(&fri_transcript);
            let mut beta_bytes = [0u8; 32];
            beta_bytes.copy_from_slice(&fri_hash[..32]);
            let beta = Fr::from_le_bytes_mod_order(&beta_bytes);
            d_evals = fri_fold(&d_evals, &xs, beta);
            xs.truncate(d_evals.len()); // keep only the first half
            for x in &mut xs {
                *x = *x * *x;
            }
            if d_evals.iter().all(|v| *v == d_evals[0]) {
                break;
            }

            folding_steps += 1;
        }

        println!("Constraint degree: {}", &c_poly.degree());
        println!("Quotient degree: {}", &q_poly.degree());
        println!("DEEP degree: {}", &d_poly_degree);
        println!("Folding steps: {}", &folding_steps);

        assert_eq!(folding_steps, 5);
        StarkProof { fri_challenges }
    }
}

fn get_random_z(
    extended_domain: &GeneralEvaluationDomain<Fr>,
    shifted_domain: &GeneralEvaluationDomain<Fr>,
) -> Fr {
    let mut rng = thread_rng();
    let ext_set: std::collections::HashSet<Fr> = extended_domain.elements().collect();
    let shift_set: std::collections::HashSet<Fr> = shifted_domain.elements().collect();
    let g = extended_domain.group_gen();

    loop {
        let z = Fr::rand(&mut rng);
        if !ext_set.contains(&z)
            && !shift_set.contains(&z)
            && !shift_set.contains(&(g * z))
            && !shift_set.contains(&(g * g * z))
        {
            return z;
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{program::trace::ExecutionTrace, prover::StarkProver};
    use ark_bls12_381::Fr;

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
        for i in 1..50 {
            trace[i] = i as u64 * 3143;
        }
        let trace_field: Vec<Fr> = trace.iter().map(|x| Fr::from(*x)).collect();
        execution_trace.insert_column(trace_field);
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
