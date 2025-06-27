#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::Field;
    use toyni::{prover::StarkProver, verifier::StarkVerifier, program::{constraints::ConstraintSystem, trace::ExecutionTrace}};
    use std::collections::HashMap;
    use toyni::prover::build_proof_transcript;
    use ark_poly::domain::GeneralEvaluationDomain;
    use ark_poly::EvaluationDomain;
    use toyni::prover::generate_spot_check_challenges;

    #[test]
    fn test_valid_proof() {
        let mut trace = ExecutionTrace::new(4, 1);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(trace.clone(), constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_invalid_proof() {
        let mut trace = ExecutionTrace::new(4, 1);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i + 1); // invalid
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(trace.clone(), constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(trace.height as usize);
        assert!(!verifier.verify(&proof));
    }

    #[test]
    fn test_larger_trace() {
        let mut trace = ExecutionTrace::new(8, 1);
        for i in 0..8 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(trace.clone(), constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_multiple_variables() {
        let mut trace = ExecutionTrace::new(4, 2);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            row.insert("y".to_string(), i * 2);
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        // x[n+1] = x[n] + 1
        constraints.add_transition_constraint(
            "increment_x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        // y[n] = 2 * x[n]
        constraints.add_transition_constraint(
            "y_is_double_x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, _| {
                let x = Fr::from(*current.get("x").unwrap());
                let y = Fr::from(*current.get("y").unwrap());
                y - x * Fr::from(2u64)
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(trace.clone(), constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_zero_values() {
        let mut trace = ExecutionTrace::new(4, 1);
        for _ in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), 0); // All zeros
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "zero_sequence".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n // Should be zero
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_zero".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(trace.clone(), constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_complex_constraints() {
        let mut trace = ExecutionTrace::new(4, 2);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            row.insert("y".to_string(), i * i); // y = x^2
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        // x[n+1] = x[n] + 1
        constraints.add_transition_constraint(
            "increment_x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        // y[n] = x[n]^2
        constraints.add_transition_constraint(
            "y_is_x_squared".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, _| {
                let x = Fr::from(*current.get("x").unwrap());
                let y = Fr::from(*current.get("y").unwrap());
                y - x * x
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(trace.clone(), constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(trace.height as usize);
        assert!(verifier.verify(&proof));
    }

    #[test]
    fn test_invalid_complex_constraints() {
        let mut trace = ExecutionTrace::new(4, 2);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            row.insert("y".to_string(), i * i + 1); // y = x^2 + 1 (invalid)
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        // x[n+1] = x[n] + 1
        constraints.add_transition_constraint(
            "increment_x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        // y[n] = x[n]^2
        constraints.add_transition_constraint(
            "y_is_x_squared".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, _| {
                let x = Fr::from(*current.get("x").unwrap());
                let y = Fr::from(*current.get("y").unwrap());
                y - x * x
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(trace.clone(), constraints);
        let proof = prover.generate_proof();
        let verifier = StarkVerifier::new(trace.height as usize);
        assert!(!verifier.verify(&proof));
    }

    #[test]
    fn test_challenge_consistency() {
        let mut trace = ExecutionTrace::new(4, 1);
        for i in 0..4 {
            let mut row = HashMap::new();
            row.insert("x".to_string(), i);
            trace.insert_column(row);
        }

        let mut constraints = ConstraintSystem::default();
        constraints.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_n = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - x_n - Fr::ONE
            }),
        );
        constraints.add_boundary_constraint(
            "starts_at_0".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let prover = StarkProver::new(trace.clone(), constraints);
        let proof = prover.generate_proof();
        println!("Prover FRI challenges: {:?}", proof.fri_challenges);
        println!("Prover spot-check challenges: {:?}", proof.verifier_random_challenges);

        // Reconstruct the verifier's challenges using the same transcript logic
        let proof_transcript = build_proof_transcript(
            proof.fri_layers.last().unwrap(),
            &proof.fri_layers,
            &proof.fri_challenges,
            &proof.combined_constraint,
            &proof.folding_commitment_trees,
        );
        let extended_domain = GeneralEvaluationDomain::<Fr>::new(trace.height as usize * 8).unwrap();
        let verifier_random_challenges = generate_spot_check_challenges(
            &proof_transcript,
            &extended_domain,
            80,
        );
        assert_eq!(proof.verifier_random_challenges, verifier_random_challenges);
    }
}
