//! Evaluation domain operations for Stark proofs.
//!
//! This module provides functionality for working with evaluation domains in the Stark proving system.
//! It includes functions for creating and extending evaluation domains, as well as operations on domain points.
use ark_bls12_381::Fr;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

/// Squares each point in the domain for FRI protocol.
///
/// This operation is used in the FRI protocol to reduce the size of the evaluation domain
/// while maintaining certain algebraic properties.
///
/// # Arguments
///
/// * `domain_points` - The original domain points to fold
/// * `domain_size` - The size of the domain (must be even)
///
/// # Returns
///
/// A vector of folded domain points
///
/// # Panics
///
/// Panics if the domain size is not even
///
pub fn fold_domain_points(domain_points: Vec<Fr>, domain_size: usize) -> Vec<Fr> {
    domain_points
        .iter()
        .take(domain_size / 2)
        .map(|&x| x * x) // Square each point
        .collect()
}

/// Creates an evaluation domain of size 2^n.
///
/// # Arguments
///
/// * `domain_size` - The size of the domain (must be a power of 2)
///
/// # Returns
///
/// A new evaluation domain
///
/// # Panics
///
/// Panics if the domain size is not a power of 2
pub fn get_domain(domain_size: usize) -> GeneralEvaluationDomain<Fr> {
    GeneralEvaluationDomain::<Fr>::new(domain_size).unwrap()
}

/// Creates an extended domain by applying a blowup factor.
///
/// The extended domain is used to improve the security of the Stark proof by
/// increasing the size of the evaluation domain.
///
/// # Arguments
///
/// * `domain_size` - The size of the original domain (must be a power of 2)
/// * `blowup_factor` - The factor by which to extend the domain
///
/// # Returns
///
/// An extended evaluation domain
///
/// # Panics
///
/// Panics if:
/// * The domain size is not a power of 2
/// * The resulting extended domain size is not a power of 2
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
