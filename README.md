# toyni

A modular STARK library with **two proving engines behind one generic
frontend**: the classical roots-of-unity STARK over BabyBear, and a circle
STARK over Mersenne-31. The same AIR definition, example programs, Merkle
commitments and Fiat–Shamir transcript run on either engine unchanged.

> [!CAUTION]
> Research code. Not audited; do not use where a broken proof would have
> real-world consequences.

## Two engines, one architecture

```
                    ┌───────────────────────────────────────────┐
   your program ───►│  Air<B>  (transition + boundary, degree ≤2)│
                    │  TraceTable<B>   — generic over base field │
                    └───────────────────┬───────────────────────┘
                                        │  StarkProver / StarkVerifier
                          ┌─────────────┴──────────────┐
                          ▼                            ▼
              ┌───────────────────────┐   ┌──────────────────────────┐
   Engine ──► │  ClassicalEngine      │   │  CircleEngine  (feature)  │
   trait      │  BabyBear (2^31−2^27+1)│   │  M31 (2^31−1)             │
              │  NTT · coset domains  │   │  circle FFT · twin-cosets │
              │  2-adic FRI           │   │  circle FRI (λ-decomp)    │
              └───────────────────────┘   │  optional Metal GPU       │
                          │                └──────────────────────────┘
                          └────────► shared: field traits, Merkle (SHA-256),
                                     transcript, DEEP-ALI prover/verifier
```

Everything that differs between roots-of-unity and circle geometry — domain
points, FFT, vanishing polynomials, single-point quotient denominators, the
out-of-domain point type, and the FRI fold — lives behind the [`Engine`]
trait (`src/engine/mod.rs`). Both engines realize the same DEEP-ALI protocol:
commit the trace, commit the constraint composition, sample out-of-domain,
batch single-point DEEP quotients, FRI to a constant final layer.

```rust
use toyni::{examples::FibonacciAir, ClassicalEngine, ProofOptions,
            StarkProver, StarkVerifier};
use toyni::field::BabyBear;

let (air, trace) = FibonacciAir::<BabyBear>::with_trace(16);
let o = ProofOptions::default();
let proof = StarkProver::<ClassicalEngine, _>::new(&air, o).prove(&trace);
assert!(StarkVerifier::<ClassicalEngine, _>::new(&air, o).verify(&proof));
```

Swap `ClassicalEngine` → `CircleEngine` and `BabyBear` → `M31` to run the same
AIR on the circle engine (with `--features circle`). The example AIRs
(`FibonacciAir`, `RangeCheckAir`, `HashChainAir`) are generic over the base
field and run on both.

## Layout

```
src/
  field/        Field / ExtField / ConstraintField traits + both towers:
                babybear, babybear_ext | m31, cm31, qm31  (circle feature)
  air.rs        Air trait, TraceTable, Boundary  (the generic frontend)
  merkle.rs     tagged SHA-256 Merkle (shared)
  transcript.rs Fiat–Shamir byte sponge (shared)
  proof.rs      ProofOptions, StarkProof<E>
  prover.rs     generic DEEP-ALI prover over any Engine
  verifier.rs   generic verifier
  engine/
    mod.rs      the Engine trait + Mode
    classical/  ClassicalEngine: ntt, coset domains, 2-adic FRI
    circle/     CircleEngine: geometry, cfft, circle FRI, backend, metal
  examples/     fibonacci, range_check, hash_chain  (generic over B)
```

## Features

| feature | effect |
|---------|--------|
| `parallel` (default) | rayon-parallel FFT / Merkle / composition |
| `circle` | compile the circle engine + M31 field tower |
| `metal` | circle engine's Metal GPU backend (implies `circle`) |
| `cuda` | classical NTT on CUDA (existing) |

## Benchmarks: roots-of-unity vs circle

Same AIR, same protocol parameters (rate `2^-3`, 44 queries, zero-knowledge),
same generic prover — only the engine differs. Apple M3 Max, multicore CPU:

| workload | roots-of-unity prove | circle prove | speed-up | verify (RoU / circle) |
|----------|---------------------:|-------------:|---------:|----------------------:|
| Fibonacci 2^12       | 84 ms  | 75 ms  | 1.1× | 5.1 / 4.1 ms |
| Fibonacci 2^14       | 318 ms | 266 ms | 1.2× | 7.0 / 5.4 ms |
| Fibonacci 2^16       | 1.28 s | 1.04 s | 1.2× | 12 / 7.3 ms |
| hash-chain 2^12 × 24 | 510 ms | 255 ms | **2.0×** | 5.2 / 4.4 ms |
| hash-chain 2^14 × 24 | 2.17 s | 965 ms | **2.2×** | 7.1 / 5.5 ms |

Circle wins on the wide multi-column workload where M31's faster arithmetic
(shift-and-add reduction vs BabyBear's multiply-reduce) and the halved FRI
blowup dominate. The circle engine additionally has a Metal GPU backend
(`--features metal`) that is bit-identical to the CPU path and gives a further
~3× on large traces (see `src/engine/circle/metal`).

Reproduce: `cargo run --release --features circle --example bench`.

## Tests

`cargo test --features circle` — 67 library unit tests + 12 end-to-end tests
that prove/verify every example AIR on **both** engines, including tamper
rejection (mutating any commitment, opening, OOD value, or FRI layer is
rejected). Clippy clean on every feature combination.

## Used by

[toyni-zkvm](https://github.com/jonas089/toyni-zkvm) builds a field-native
RISC-style zkVM on toyni's primitives (LogUp memory/register consistency,
grand-product permutation arguments). It runs on the roots-of-unity engine;
the unified engine architecture here is what makes a circle-engine zkVM a
drop-in change of the geometry layer rather than a rewrite.
