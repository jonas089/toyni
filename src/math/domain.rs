use ark_bls12_381::Fr;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

pub fn fold_domain_points(domain_points: Vec<Fr>, domain_size: usize) -> Vec<Fr> {
    domain_points
        .iter()
        .take(domain_size / 2)
        .map(|&x| x * x) // Square each point
        .collect()
}

pub fn get_domain(domain_size: usize) -> GeneralEvaluationDomain<Fr> {
    GeneralEvaluationDomain::<Fr>::new(domain_size).unwrap()
}

pub fn get_extended_domain(
    domain_size: usize,
    blowup_factor: usize,
) -> GeneralEvaluationDomain<Fr> {
    GeneralEvaluationDomain::<Fr>::new(domain_size * blowup_factor).unwrap()
}

#[test]
fn test_general_evaluation_domain() {
    let original_domain_size = 4;
    let blowup_factor = 8;
    let domain = get_domain(original_domain_size);
    let extended_domain = get_extended_domain(original_domain_size, blowup_factor);

    // Verify domain sizes
    assert_eq!(domain.size(), 4);
    assert_eq!(extended_domain.size(), 4 * blowup_factor);

    // Verify that the extended domain contains the original domain points
    let original_points: Vec<Fr> = domain.elements().collect();
    let extended_points: Vec<Fr> = extended_domain.elements().collect();

    // The original points should appear in the extended domain at regular intervals
    for (i, &point) in original_points.iter().enumerate() {
        assert_eq!(point, extended_points[i * blowup_factor]);
    }
}
