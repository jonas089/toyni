use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::merkle::MerkleTree;
use crate::{digest_sha2, program::trace::ExecutionTrace};
use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, BigInteger, Field, PrimeField, UniformRand};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial};
use rand::{Rng, thread_rng};

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
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace_len * 8).unwrap();

        let z_poly = ToyniPolynomial::from_dense_poly(domain.vanishing_polynomial().into());
        let r_poly = random_poly(2);

        fn fibonacci_constraint(t2: Fr, t1: Fr, t0: Fr) -> Fr {
            t2 - (t1 + t0)
        }

        let trace_poly = self
            .trace
            .interpolate_column(&domain.elements().collect::<Vec<Fr>>(), 0);

        let g = domain.group_gen();

        let extended_and_shifted_domain: Vec<Fr> = extended_domain
            .elements()
            .map(|x| x * Fr::from(7))
            .collect();

        let z = get_random_z(&domain, &extended_domain, Fr::from(7));
        let tz_0 = trace_poly.evaluate(z);
        let tz_1 = trace_poly.evaluate(g * z);
        let tz_2 = trace_poly.evaluate(g * g * z);

        let mut d_evals = vec![];
        for element in extended_and_shifted_domain {
            let c_z = fibonacci_constraint(tz_2, tz_1, tz_0);
            let element = element;
            let t0 = trace_poly.evaluate(element);
            let t1 = trace_poly.evaluate(g * element);
            let t2 = trace_poly.evaluate(g * g * element);
            let c_x = fibonacci_constraint(t2, t1, t0);

            let d_x =
                (c_x - c_z) / (element - z) + r_poly.evaluate(element) * z_poly.evaluate(element);

            if d_x == Fr::ZERO {
                continue;
            }
            d_evals.push(d_x);
        }

        let mut fri_layers = vec![d_evals.clone()];
        let mut fri_challenges = Vec::new();
        let mut folding_commitment_trees: Vec<MerkleTree> = Vec::new();
        let mut fri_transcript = Vec::new();

        // folding the D(x) polynomial
        let mut folding_steps = 0;
        let mut fri_final_constant = false;

        while d_evals.len() > 1 {
            folding_steps += 1;

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

fn get_random_z(
    domain: &GeneralEvaluationDomain<Fr>,
    extended_domain: &GeneralEvaluationDomain<Fr>,
    shift: Fr,
) -> Fr {
    let mut rng = thread_rng();
    let extended_shifted: Vec<Fr> = extended_domain.elements().map(|x| x * shift).collect();

    loop {
        let z = Fr::rand(&mut rng);
        if !extended_shifted.contains(&z)
            && !extended_shifted.contains(&(domain.group_gen() * z))
            && !extended_shifted.contains(&(domain.group_gen().square() * z))
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
        let trace: Vec<i64> = vec![
            1,
            1,
            2,
            3,
            5,
            8,
            13,
            21,
            34,
            56,
            89,
            144,
            233,
            377,
            610,
            987,
            1597,
            2584,
            4181,
            6765,
            10946,
            17711,
            28657,
            46368,
            75025,
            121393,
            196418,
            317811,
            514229,
            832040,
            1346269,
            2178309,
            3524578,
            5702887,
            9227465,
            14930352,
            24157817,
            39088169,
            63245986,
            102334155,
            165580141,
            267914296,
            433494437,
            701408733,
            1134903170,
            1836311903,
            2971215073,
            4807526976,
            7778742049,
            12586269025,
            20365011074,
            32951280099,
            53316291173,
            86267571272,
            139583862445,
            225851433717,
            365435296162,
            591286729879,
            956722026041,
            1548008755920,
            2504730781961,
            4052739537881,
            6557470319842,
            10610209857723,
        ];
        let trace_field: Vec<Fr> = trace.iter().map(|x| Fr::from(*x)).collect();
        execution_trace.insert_column(trace_field);
        let stark = StarkProver::new(execution_trace.clone());
        let proof = stark.generate_proof();
        /*let verifier = StarkVerifier::new(execution_trace.trace.len());
        assert!(verifier.verify(&proof));*/
    }

    #[test]
    fn test_invalid_fibonacci() {
        let mut execution_trace = ExecutionTrace::new();
        let trace: Vec<i64> = vec![
            34563456, 354634, 546345, 5346, 3456, 5436453, 3456345, 3546345, 34, 56, 1344, 144,
            233, 377, 610, 555, 1597, 2584, 4181, 6666, 10946, 17711, 45345, 46368, 75025, 121393,
            2462432, 317811, 514229, 832040, 346234, 2178309, 3524578, 5702887, 24352345, 14930352,
            24157817, 34564356, 234523, 64536435, 234523, 34434, 4453, 7536, 3546, 34567, 3456,
            4563456, 56456, 3745, 3456, 88456, 3423, 6677, 555555, 234234, 23444, 234, 234,
            2314123, 3434, 5555555, 34234, 123,
        ];
        let trace_field: Vec<Fr> = trace.iter().map(|x| Fr::from(*x)).collect();
        execution_trace.insert_column(trace_field);
        let stark = StarkProver::new(execution_trace.clone());
        let proof = stark.generate_proof();
        /*let verifier = StarkVerifier::new(execution_trace.trace.len());
        assert!(verifier.verify(&proof));*/
    }

    #[test]
    fn fibonacci_list() {
        let n = 64;
        let mut fibs: Vec<u64> = Vec::with_capacity(n);
        let mut a = 1;
        let mut b = 1;

        for _ in 0..n {
            fibs.push(a);
            let next = a + b;
            a = b;
            b = next;
        }

        println!("fibs: {:?}", &fibs);
    }
}
