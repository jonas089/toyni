use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_poly::{
    EvaluationDomain, Evaluations, GeneralEvaluationDomain, univariate::DensePolynomial,
};

pub fn fri_fold(evals: &[Fr], xs: &[Fr], beta: Fr) -> Vec<Fr> {
    assert!(evals.len() % 2 == 0, "Evaluations length must be even");
    let half = evals.len() / 2;
    let half_inv = Fr::from(2u64).inverse().unwrap();
    let mut result = Vec::with_capacity(half);

    for i in 0..half {
        let a = evals[i]; // f(xᵢ)
        let b = evals[i + half]; // f(-xᵢ)
        let x = xs[i]; // xᵢ

        // average part
        let avg = (a + b) * half_inv;
        // “odd” part, scaled and then divided by x to push onto x²‐domain
        let diff = (a - b) * half_inv;
        let folded = avg + diff * beta * x.inverse().unwrap();

        result.push(folded);
    }

    result
}

pub fn interpolate_poly(xs: &[Fr], ys: &[Fr]) -> DensePolynomial<Fr> {
    assert_eq!(xs.len(), ys.len(), "Mismatched lengths");
    let domain =
        GeneralEvaluationDomain::<Fr>::new(xs.len()).expect("Domain size must be a power of 2");
    let evals = Evaluations::from_vec_and_domain(ys.to_vec(), domain);
    evals.interpolate()
}
