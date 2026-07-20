//! Benchmark both engines on the same AIRs. Run:
//!   cargo run --release --features circle --example bench
//!   cargo run --release --features metal  --example bench

use std::time::Instant;

use toyni::examples::{FibonacciAir, HashChainAir};
use toyni::field::{BabyBear, Field};
use toyni::{Air, ClassicalEngine, Engine, ProofOptions, StarkProof, StarkProver, StarkVerifier, TraceTable};

fn run<E: Engine, A: Air<E::Base>>(name: &str, air: &A, trace: &TraceTable<E::Base>) {
    let o = ProofOptions::default();
    let prover = StarkProver::<E, _>::new(air, o);
    let _ = prover.prove(trace); // warmup (shader compile, etc.)
    let t0 = Instant::now();
    let proof: StarkProof<E> = prover.prove(trace);
    let dt = t0.elapsed();
    let t1 = Instant::now();
    let ok = StarkVerifier::<E, _>::new(air, o).verify(&proof);
    println!(
        "  {name:<26} [{}]  prove {dt:>10.2?}  verify {:>8.2?}  (ok={ok})",
        E::label(),
        t1.elapsed()
    );
}

fn main() {
    let backend = if toyni::metal_available() { "metal" } else { "cpu" };
    println!("toyni engine benchmark  (circle GPU backend: {backend})\n");

    for log in [12u32, 14, 16] {
        println!("Fibonacci 2^{log}:");
        let (air_c, trace_c) = FibonacciAir::<BabyBear>::with_trace(log);
        run::<ClassicalEngine, _>("fibonacci", &air_c, &trace_c);
        #[cfg(feature = "circle")]
        {
            use toyni::field::M31;
            use toyni::CircleEngine;
            let (air, trace) = FibonacciAir::<M31>::with_trace(log);
            run::<CircleEngine, _>("fibonacci", &air, &trace);
        }
    }

    for log in [12u32, 14] {
        println!("Hash-chain 2^{log} x 24 columns:");
        let (air_c, trace_c) = HashChainAir::<BabyBear>::with_trace(log, 7);
        run::<ClassicalEngine, _>("hash-chain", &air_c, &trace_c);
        #[cfg(feature = "circle")]
        {
            use toyni::field::M31;
            use toyni::CircleEngine;
            let (air, trace) = HashChainAir::<M31>::with_trace(log, 7);
            run::<CircleEngine, _>("hash-chain", &air, &trace);
        }
    }
    let _ = BabyBear::ONE;
}
