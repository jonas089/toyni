//! TOYNI-03: the Fibonacci AIR has no boundary constraints (see REPORT.md).

use toyni::babybear::BabyBear;
use toyni::fibonacci::StarkProver;
use toyni::program::trace::ExecutionTrace;
use toyni::verifier::StarkVerifier;

fn prove(seq: Vec<BabyBear>) -> bool {
    let mut t = ExecutionTrace::new();
    t.insert_column(seq);
    StarkVerifier.verify(&StarkProver::new(t).generate_proof(false))
}

fn recurrence(a: u64, b: u64, n: usize) -> Vec<BabyBear> {
    let mut v = vec![BabyBear::new(a), BabyBear::new(b)];
    while v.len() < n {
        let k = v.len();
        v.push(v[k - 1] + v[k - 2]);
    }
    v
}

#[test]
fn all_zero_sequence_verifies() {
    assert!(prove(vec![BabyBear::zero(); 64]));
}

#[test]
fn arbitrary_seeds_verify() {
    assert!(prove(recurrence(1_000_000, 7, 64)));
}

