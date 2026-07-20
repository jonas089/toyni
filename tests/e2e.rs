//! End-to-end prove/verify over both engines, plus tamper-rejection checks.
//! The same generic AIRs run on the classical (BabyBear) and circle (M31)
//! engines, exercising the unified architecture.

use toyni::examples::{FibonacciAir, HashChainAir, RangeCheckAir};
use toyni::field::{BabyBear, BabyBearExt, Field};
use toyni::{ClassicalEngine, ProofOptions, StarkProver, StarkVerifier};

fn opts() -> ProofOptions {
    ProofOptions {
        log_blowup: 2,
        num_queries: 16,
        zk: true,
    }
}

// ── classical (roots-of-unity, BabyBear) ───────────────────────────────

#[test]
fn classical_fibonacci() {
    let (air, trace) = FibonacciAir::<BabyBear>::with_trace(8);
    let o = opts();
    let proof = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
    assert!(StarkVerifier::<ClassicalEngine, _>::new(&air, o).verify(&proof));
}

#[test]
fn classical_range_check() {
    let (air, trace) = RangeCheckAir::<BabyBear>::with_trace(8, 0xFEED);
    let o = opts();
    let proof = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
    assert!(StarkVerifier::<ClassicalEngine, _>::new(&air, o).verify(&proof));
}

#[test]
fn classical_hash_chain() {
    let (air, trace) = HashChainAir::<BabyBear>::with_trace(7, 99);
    let o = opts();
    let proof = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
    assert!(StarkVerifier::<ClassicalEngine, _>::new(&air, o).verify(&proof));
}

#[test]
fn classical_default_options() {
    let (air, trace) = FibonacciAir::<BabyBear>::with_trace(10);
    let o = ProofOptions::default();
    let proof = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
    assert!(StarkVerifier::<ClassicalEngine, _>::new(&air, o).verify(&proof));
}

#[test]
fn classical_zk_randomizes() {
    let (air, trace) = FibonacciAir::<BabyBear>::with_trace(8);
    let o = opts();
    let p1 = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
    let p2 = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
    assert!(StarkVerifier::<ClassicalEngine, _>::new(&air, o).verify(&p1));
    assert!(StarkVerifier::<ClassicalEngine, _>::new(&air, o).verify(&p2));
    assert_ne!(p1.ood_trace[0][0], p2.ood_trace[0][0]);
}

#[test]
fn classical_tamper_rejected() {
    let (air, trace) = FibonacciAir::<BabyBear>::with_trace(8);
    let o = opts();
    let proof = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
    let verify = |p: &_| StarkVerifier::<ClassicalEngine, _>::new(&air, o).verify(p);
    assert!(verify(&proof));

    let mut p = proof.clone();
    p.trace_root[0] ^= 0xff;
    assert!(!verify(&p));
    let mut p = proof.clone();
    p.composition_root[0] ^= 0xff;
    assert!(!verify(&p));
    let mut p = proof.clone();
    p.ood_trace[0][0] += BabyBearExt::ONE;
    assert!(!verify(&p));
    let mut p = proof.clone();
    p.ood_composition += BabyBearExt::ONE;
    assert!(!verify(&p));
    let mut p = proof.clone();
    p.fri.final_layer[0] += BabyBearExt::ONE;
    assert!(!verify(&p));
    let mut p = proof.clone();
    p.queries[0].trace.values[0] += BabyBear::ONE;
    assert!(!verify(&p));
    let mut p = proof.clone();
    p.queries.pop();
    assert!(!verify(&p));
}

#[test]
fn classical_wrong_output_rejected() {
    let (air, trace) = FibonacciAir::<BabyBear>::with_trace(8);
    let o = opts();
    let proof = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
    let wrong = FibonacciAir::<BabyBear> {
        claimed_output: air.claimed_output + BabyBear::ONE,
    };
    assert!(!StarkVerifier::<ClassicalEngine, _>::new(&wrong, o).verify(&proof));
}

// ── circle (M31), same AIRs ────────────────────────────────────────────

#[cfg(feature = "circle")]
mod circle {
    use super::*;
    use toyni::field::{M31, QM31};
    use toyni::CircleEngine;

    #[test]
    fn circle_fibonacci() {
        let (air, trace) = FibonacciAir::<M31>::with_trace(8);
        let o = opts();
        let proof = StarkProver::<CircleEngine, _>::new(&air, o).prove(&trace);
        assert!(StarkVerifier::<CircleEngine, _>::new(&air, o).verify(&proof));
    }

    #[test]
    fn circle_range_check() {
        let (air, trace) = RangeCheckAir::<M31>::with_trace(8, 0xFEED);
        let o = opts();
        let proof = StarkProver::<CircleEngine, _>::new(&air, o).prove(&trace);
        assert!(StarkVerifier::<CircleEngine, _>::new(&air, o).verify(&proof));
    }

    #[test]
    fn circle_hash_chain() {
        let (air, trace) = HashChainAir::<M31>::with_trace(7, 99);
        let o = opts();
        let proof = StarkProver::<CircleEngine, _>::new(&air, o).prove(&trace);
        assert!(StarkVerifier::<CircleEngine, _>::new(&air, o).verify(&proof));
    }

    #[test]
    fn circle_tamper_rejected() {
        let (air, trace) = FibonacciAir::<M31>::with_trace(8);
        let o = opts();
        let proof = StarkProver::<CircleEngine, _>::new(&air, o).prove(&trace);
        let verify = |p: &_| StarkVerifier::<CircleEngine, _>::new(&air, o).verify(p);
        assert!(verify(&proof));

        let mut p = proof.clone();
        p.trace_root[0] ^= 0xff;
        assert!(!verify(&p));
        let mut p = proof.clone();
        p.ood_composition += QM31::ONE;
        assert!(!verify(&p));
        let mut p = proof.clone();
        p.fri.lambda += QM31::ONE;
        assert!(!verify(&p));
        let mut p = proof.clone();
        p.queries[0].trace.values[0] += M31::ONE;
        assert!(!verify(&p));
    }

    #[test]
    fn circle_zk_randomizes() {
        let (air, trace) = FibonacciAir::<M31>::with_trace(8);
        let o = opts();
        let p1 = StarkProver::<CircleEngine, _>::new(&air, o).prove(&trace);
        let p2 = StarkProver::<CircleEngine, _>::new(&air, o).prove(&trace);
        assert!(StarkVerifier::<CircleEngine, _>::new(&air, o).verify(&p1));
        assert_ne!(p1.ood_trace[0][0], p2.ood_trace[0][0]);
    }
}
