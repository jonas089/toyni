//! Fast Reed-Solomon Interactive Oracle Proof (FRI) implementation.
//!
//! The FRI protocol proves polynomial low-degree by iteratively folding evaluations
//! and reducing domain size. Each round combines pairs of points using a random challenge.

use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_poly::{
    EvaluationDomain, Evaluations, GeneralEvaluationDomain, univariate::DensePolynomial,
};

/// Folds evaluations using FRI protocol with challenge beta.
pub fn fri_fold(evals: &[Fr], beta: Fr) -> Vec<Fr> {
    assert!(evals.len() % 2 == 0, "Evaluations length must be even");
    let mut result = Vec::with_capacity(evals.len() / 2);
    let half = evals.len() / 2;
    let half_inv = Fr::from(2u64).inverse().unwrap();

    for i in 0..half {
        let a = evals[i];
        let b = evals[i + half];
        let folded = (a + b) * half_inv + (a - b) * half_inv * beta;
        result.push(folded);
    }

    result
}

/// Interpolates polynomial from points using FFT.
pub fn interpolate_poly(xs: &[Fr], ys: &[Fr]) -> DensePolynomial<Fr> {
    assert_eq!(xs.len(), ys.len(), "Mismatched lengths");
    let domain =
        GeneralEvaluationDomain::<Fr>::new(xs.len()).expect("Domain size must be a power of 2");
    let evals = Evaluations::from_vec_and_domain(ys.to_vec(), domain);
    evals.interpolate()
}
