#[cfg(test)]
mod tests {
    use toyni::babybear::BabyBear;
    use toyni::math::{
        domain::BabyBearDomain,
        fri::fri_fold,
        polynomial::Polynomial,
    };

    #[test]
    fn test_domain_properties() {
        let blowup_factor = 8;
        let domain = BabyBearDomain::new(4);
        let extended_domain = BabyBearDomain::new(4 * blowup_factor);

        assert_eq!(domain.size(), 4);
        assert_eq!(extended_domain.size(), 4 * blowup_factor);

        let original_points = domain.elements();
        let extended_points = extended_domain.elements();

        for (i, &point) in original_points.iter().enumerate() {
            assert_eq!(point.value, extended_points[i * blowup_factor].value);
        }
    }

    #[test]
    fn test_polynomial_division() {
        // Test case: (x^3 + 2x^2 + 3x + 4) / (x + 1)
        let dividend = Polynomial::new(vec![
            BabyBear::new(4),
            BabyBear::new(3),
            BabyBear::new(2),
            BabyBear::new(1),
        ]);
        let divisor = Polynomial::new(vec![BabyBear::new(1), BabyBear::new(1)]);
        let result = dividend.divide(&divisor).unwrap();
        let (quotient, remainder) = result;
        assert_eq!(quotient.coefficients[0].value, 2);
        assert_eq!(quotient.coefficients[1].value, 1);
        assert_eq!(quotient.coefficients[2].value, 1);
        assert_eq!(remainder.coefficients[0].value, 2);
    }

    #[test]
    fn test_polynomial_division_exact() {
        // Test case: (x^2 + 2x + 1) / (x + 1)
        let dividend = Polynomial::new(vec![
            BabyBear::new(1),
            BabyBear::new(2),
            BabyBear::new(1),
        ]);
        let divisor = Polynomial::new(vec![BabyBear::new(1), BabyBear::new(1)]);
        let result = dividend.divide(&divisor).unwrap();
        let (quotient, remainder) = result;
        assert_eq!(quotient.coefficients[0].value, 1);
        assert_eq!(quotient.coefficients[1].value, 1);
        assert!(remainder.coefficients.is_empty() || remainder.coefficients[0].is_zero());
    }

    #[test]
    fn test_polynomial_division_zero() {
        let dividend = Polynomial::new(vec![
            BabyBear::new(1),
            BabyBear::new(2),
            BabyBear::new(1),
        ]);
        let divisor = Polynomial::new(vec![BabyBear::zero()]);
        assert!(dividend.divide(&divisor).is_none());
    }

    #[test]
    fn test_polynomial_addition() {
        let p1 = Polynomial::new(vec![BabyBear::new(1), BabyBear::new(2), BabyBear::new(3)]);
        let p2 = Polynomial::new(vec![BabyBear::new(4), BabyBear::new(5), BabyBear::new(6)]);
        let sum = p1.add(&p2);
        assert_eq!(sum.coefficients[0].value, 5);
        assert_eq!(sum.coefficients[1].value, 7);
        assert_eq!(sum.coefficients[2].value, 9);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = Polynomial::new(vec![BabyBear::new(1), BabyBear::new(2)]); // 2x + 1
        let p2 = Polynomial::new(vec![BabyBear::new(3), BabyBear::new(4)]); // 4x + 3
        let product = p1.multiply(&p2);
        assert_eq!(product.coefficients[0].value, 3);
        assert_eq!(product.coefficients[1].value, 10);
        assert_eq!(product.coefficients[2].value, 8);
    }

    #[test]
    fn test_polynomial_multiplication_zero() {
        let p1 = Polynomial::new(vec![BabyBear::new(1), BabyBear::new(2)]);
        let p2 = Polynomial::new(vec![BabyBear::zero()]);
        let product = p1.multiply(&p2);
        assert!(product.coefficients.is_empty());
    }

    #[test]
    fn test_fri_folding() {
        let domain = BabyBearDomain::new(4);

        // Create polynomial: 1 + 2x + x^2
        let coeffs = vec![BabyBear::new(1), BabyBear::new(2), BabyBear::new(1)];

        // Evaluate on domain
        let eval_vec = domain.fft(&coeffs);

        // Random beta
        let mut rng = rand::thread_rng();
        let beta = BabyBear::random(&mut rng);

        // Fold
        let xs = domain.elements();
        let folded_evals = fri_fold(&eval_vec, &xs, beta);

        // Folded domain (squared points)
        let folded_domain_points: Vec<BabyBear> = xs
            .iter()
            .take(2)
            .map(|&x| x * x)
            .collect();

        // Verify folded evaluations match an interpolated polynomial
        let folded_poly = Polynomial::lagrange_interpolate(&folded_domain_points, &folded_evals);

        for (i, &eval) in folded_evals.iter().enumerate() {
            let x = folded_domain_points[i];
            let y = folded_poly.evaluate(x);
            assert_eq!(eval.value, y.value);
        }
    }
}
