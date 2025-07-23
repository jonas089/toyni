#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use toyni::{program::trace::ExecutionTrace, prover::StarkProver};

    #[test]
    fn test_constraint_poly() {
        let mut execution_trace = ExecutionTrace::new();
        execution_trace.insert_column(vec![
            Fr::from(2),
            Fr::from(4),
            Fr::from(16),
            Fr::from(16 * 16),
        ]);
        let stark = StarkProver::new(execution_trace);
        let proof = stark.generate_proof();
        println!("{:?}", proof);
    }
}
