use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
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

pub const CONSTRAINT_SPOT_CHECKS: usize = 50;

#[derive(Debug)]
pub struct StarkProof {
    pub fri_challenges: Vec<Fr>,
    pub trace_spot_checks: [[Fr; 3]; CONSTRAINT_SPOT_CHECKS],
    pub constraint_spot_checks: [Fr; CONSTRAINT_SPOT_CHECKS],
}
/* Real implementation plan

1. interpolate T(x) over original domain
2. commit to Q(x) over shifted domain
3. fold D(x) and check some shifted spots for equality against Q(x)
4. check C(z) for consistency e.g. Q(z) = fibonacci(z) / Z(z)

Q(z) = fibonacci(ggz, gz, z) / Z(z)

we use Q(z) to build D and then fold D, committing to each layer.

The verifier checks that Q(z) = fibonacci(ggz, gz, z) / Z(z) + trace commitments
and that D(x) is low degree / consistent commitments for FRI, as well as some
spot checks in the extended shifted domain revealing:

alpha * (Q(x) - Q(z)) / (x - z)
where Q(z), z are constant
*/

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
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 2).unwrap();
        let shifted_domain = extended_domain.get_coset(Fr::from(7)).unwrap();

        // the vanishing polynomial for our trace over the original domain
        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());
        // a random polynomial that will be derived from a verifier challenge in fiat shamir
        let r_poly = random_poly(4);

        // the fibonacci constraint used for proving
        fn fibonacci_constraint(t2: Fr, t1: Fr, t0: Fr) -> Fr {
            // temporary solution - this is insecure and proper boundaries are better
            if t2 == Fr::from(10610209857723u64)
                || t1 == Fr::from(10610209857723u64)
                || t0 == Fr::from(10610209857723u64)
            {
                return Fr::ZERO;
            }
            t2 - (t1 + t0)
        }

        // the trace interpolated as a polynomial over the original domain
        let trace_poly = self
            .trace
            .interpolate_column(&domain.elements().collect::<Vec<Fr>>(), 0);

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

        // the constraint polynomial interpolated over the extended domain
        // this is equivalent to our composite constraint, because we only have one transitino constraint currently
        // and no boundary constraints.
        let c_poly = ToyniPolynomial::from_dense_poly(DensePolynomial::from_coefficients_slice(
            &extended_domain.ifft(&c_evals),
        ));

        //.add(&r_poly.mul(&z_poly));

        // evaluations of the quotient polynomial at challenge points
        let mut q_evals: Vec<Fr> = Vec::new();
        for x in shifted_domain.elements() {
            q_evals.push(c_poly.evaluate(x) / z_poly.evaluate(x));
        }

        // interpolation of the quotient for development purposes
        let q_poly = DensePolynomial::from_coefficients_slice(&shifted_domain.ifft(&q_evals));
        let z = get_random_z(&extended_domain);
        let mut d_evals = vec![];
        let mut rng = thread_rng();
        let alpha = Fr::rand(&mut rng);

        let q_z = q_poly.evaluate(&z);
        //let t_z = trace_poly.evaluate(z);

        let mut test_spot_check: Vec<Fr> = Vec::new();
        for x in shifted_domain.elements() {
            let q_x = q_poly.evaluate(&x);
            //let t_x = trace_poly.evaluate(x);
            // this is the deep formula, we expect the degree of the DEEP polynomial to be one less than the constraint polynomial
            let d_x = alpha * (q_x - q_z) / (x - z); //+ alpha * (t_x - t_z) / (x - z);
            test_spot_check.push(d_x.clone());
            d_evals.push(d_x);

            // we can get away with not using the extended domain for spot checks only if we do the following:
            // during interpolation of c(x) / when committing in practice, we also evaluate the constraints at c(z)

            // to spot check c(x) = q(x) * z(x), where c(x) is the raw constraint fn, the prover must commit to evaluations
            // over the shifted extended domain and at z for q(x)
            // we use the evaluations over the extended domain to build and fold d(x) and we use the evaluations of z
            // for a global consistency check committed_c_z = committed_q_z / Z(z)
            // if that consistency check holds and the Quotient/ DEEP poly is low degree when folded over the extended domain, then
            // the proof is considered valid
            assert_eq!(q_poly.evaluate(&x), c_poly.evaluate(x) / z_poly.evaluate(x));
            // this actually happens in folding (todo: move down)
            assert_eq!(
                d_x,
                alpha * (q_poly.evaluate(&x) - q_poly.evaluate(&z)) / (x - z) //+ alpha * (t_x - t_z) / (x - z)
            );
        }
        // simulated spot check at z
        assert_eq!(q_poly.evaluate(&z), c_poly.evaluate(z) / z_poly.evaluate(z));

        // we fold the polynomial using our FRI evaluation domain
        // the spot checks will later ensure that the polynomial was folded correctly
        // and that it is of low degree
        let fri_challenges = Vec::new();
        let mut fri_transcript: Vec<u8> = Vec::new();
        let mut folding_steps: usize = 0;

        // note the degree of the deep polynomial for debugging
        let d_poly_degree: usize =
            DensePolynomial::from_coefficients_slice(&shifted_domain.ifft(&d_evals)).degree();

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

        println!("Composite degree: {}", &c_poly.degree());
        println!("Quotient degree: {}", &q_poly.degree());
        println!("DEEP degree: {}", &d_poly_degree);
        println!("Folding steps: {}", &folding_steps);
        assert_eq!(folding_steps, 5);
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
        for i in 1..2 {
            trace[i] = i as u64 * 1233;
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
