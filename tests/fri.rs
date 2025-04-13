#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::{AdditiveGroup, Field, UniformRand};
    use ark_poly::{
        DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial,
    };
    use ark_std::test_rng;
    use toyni::math::{
        fri::{fri_fold, interpolate_poly},
        polynomial::Polynomial,
    };

    #[test]
    fn test_general_evaluation_domain() {
        let blowup_factor = 8;
        let domain = GeneralEvaluationDomain::<Fr>::new(4).unwrap();
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(4 * blowup_factor).unwrap();

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

    #[test]
    fn test_polynomial_division() {
        // Test case: (x^3 + 2x^2 + 3x + 4) / (x + 1)
        let dividend = Polynomial::new(vec![
            Fr::from(4u64),
            Fr::from(3u64),
            Fr::from(2u64),
            Fr::from(1u64),
        ]);
        let divisor = Polynomial::new(vec![Fr::from(1u64), Fr::from(1u64)]);
        let result = dividend.divide(&divisor).unwrap();
        let (quotient, remainder) = result;
        assert_eq!(quotient.coefficients[0], Fr::from(2u64));
        assert_eq!(quotient.coefficients[1], Fr::from(1u64));
        assert_eq!(quotient.coefficients[2], Fr::from(1u64));
        assert_eq!(remainder.coefficients[0], Fr::from(2u64));
    }

    #[test]
    fn test_polynomial_division_exact() {
        // Test case: (x^2 + 2x + 1) / (x + 1)
        let dividend = Polynomial::new(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(1u64)]);
        let divisor = Polynomial::new(vec![Fr::from(1u64), Fr::from(1u64)]);
        let result = dividend.divide(&divisor).unwrap();
        let (quotient, remainder) = result;
        assert_eq!(quotient.coefficients[0], Fr::from(1u64));
        assert_eq!(quotient.coefficients[1], Fr::from(1u64));
        assert!(remainder.coefficients.is_empty() || remainder.coefficients[0] == Fr::ZERO);
    }

    #[test]
    fn test_polynomial_division_zero() {
        let dividend = Polynomial::new(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(1u64)]);
        let divisor = Polynomial::new(vec![Fr::ZERO]);
        assert!(dividend.divide(&divisor).is_none());
    }

    #[test]
    fn test_polynomial_addition() {
        let p1 = Polynomial::new(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]); // 3x^2 + 2x + 1
        let p2 = Polynomial::new(vec![Fr::from(4u64), Fr::from(5u64), Fr::from(6u64)]); // 6x^2 + 5x + 4
        let sum = p1.add(&p2);
        assert_eq!(sum.coefficients[0], Fr::from(5u64));
        assert_eq!(sum.coefficients[1], Fr::from(7u64));
        assert_eq!(sum.coefficients[2], Fr::from(9u64));
    }

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = Polynomial::new(vec![Fr::from(1u64), Fr::from(2u64)]); // 2x + 1
        let p2 = Polynomial::new(vec![Fr::from(3u64), Fr::from(4u64)]); // 4x + 3
        let product = p1.multiply(&p2);
        assert_eq!(product.coefficients[0], Fr::from(3u64));
        assert_eq!(product.coefficients[1], Fr::from(10u64));
        assert_eq!(product.coefficients[2], Fr::from(8u64));
    }

    #[test]
    fn test_polynomial_multiplication_zero() {
        let p1 = Polynomial::new(vec![Fr::from(1u64), Fr::from(2u64)]); // 2x + 1
        let p2 = Polynomial::new(vec![Fr::ZERO]); // 0
        let product = p1.multiply(&p2);
        assert!(product.coefficients.is_empty());
    }

    #[test]
    fn test_fri_folding_and_interpolation() {
        // Create a test polynomial: f(x) = x^2 + 2x + 1
        let coeffs = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(1u64)];
        let poly = DensePolynomial::<Fr>::from_coefficients_vec(coeffs);
        // Create evaluation domain (size must be power of 2)
        let domain_size = 4;
        let domain = GeneralEvaluationDomain::<Fr>::new(domain_size)
            .expect("Domain size must be a power of 2");
        // Get domain points and their evaluations
        let domain_points: Vec<Fr> = domain.elements().collect();
        let eval_vec = domain.fft(&poly.coeffs);
        // Create a random beta for folding
        let mut rng = test_rng();
        let beta = Fr::rand(&mut rng);
        // Perform FRI folding on evaluations
        let folded_evals = fri_fold(&eval_vec, beta);
        // Create folded domain points by squaring the first half of the original domain
        let folded_domain_points: Vec<Fr> = domain_points
            .iter()
            .take(domain_size / 2)
            .map(|&x| x * x) // Square each point
            .collect();
        // Interpolate the folded evaluations
        let folded_poly = interpolate_poly(&folded_domain_points, &folded_evals);
        // Verify the folded polynomial has half the degree of the original
        assert_eq!(folded_poly.coeffs.len() - 1, (poly.coeffs.len() - 1) / 2);
        // Verify that the folded polynomial evaluates correctly at the folded points
        for (i, &eval) in folded_evals.iter().enumerate() {
            let x = folded_domain_points[i];
            let y = folded_poly
                .coeffs
                .iter()
                .enumerate()
                .map(|(j, &c)| c * x.pow([j as u64]))
                .sum::<Fr>();
            assert_eq!(eval, y);
        }
    }
}
