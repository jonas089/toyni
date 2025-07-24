use ark_bls12_381::Fr;
use ark_poly::{
    DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
    univariate::DensePolynomial,
};
use sha2::{Digest, Sha256};

use crate::program::trace::ExecutionTrace;

pub mod math;
pub mod merkle;
pub mod program;
pub mod prover;
pub mod verifier;

pub fn digest_sha2(data: &[u8]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher.finalize().into()
}

// experiments
pub fn experimental(trace: ExecutionTrace) {
    let domain = GeneralEvaluationDomain::<Fr>::new(trace.trace.len()).unwrap();
    let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace.trace.len() * 8).unwrap();

    fn fibonacci_constraint(ti2: Fr, ti1: Fr, ti0: Fr) -> Fr {
        ti2 - (ti1 + ti0)
    }

    let trace_poly = trace.interpolate_column(&domain.elements().collect::<Vec<Fr>>(), 0);
    let extended_points = extended_domain.elements().collect::<Vec<Fr>>();

    let trace_lde = extended_points
        .iter()
        .map(|x| trace_poly.evaluate(*x))
        .collect::<Vec<Fr>>();

    let i = 42; // choose any index < extended_points.len() - 2
    let x = extended_points[i];

    let t0 = trace_lde[i];
    let t1 = trace_lde[i + 1];
    let t2 = trace_lde[i + 2];

    let expected = fibonacci_constraint(t2, t1, t0);
    let c_poly = DensePolynomial::from_coefficients_slice(
        &extended_domain.ifft(
            &(0..(extended_points.len() - 2))
                .map(|j| {
                    let a = trace_lde[j];
                    let b = trace_lde[j + 1];
                    let c = trace_lde[j + 2];
                    fibonacci_constraint(c, b, a)
                })
                .collect::<Vec<_>>(),
        ),
    );

    let (quotient_poly, _) = c_poly.divide_by_vanishing_poly(domain);

    let actual = c_poly.evaluate(&x);

    println!("Quotient: {}", quotient_poly.degree());
    assert_eq!(
        expected, actual,
        "Constraint poly does not match direct evaluation"
    );
}

#[test]
fn experiment() {
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
    experimental(execution_trace);
}
