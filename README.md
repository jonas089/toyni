# Toyni: A STARK Implementation in Progress

> [!CAUTION]
> **This is research / hobby code.** It has **not** been audited and is not
> suitable for production use. Do not rely on it for any setting where a
> broken proof would have real-world consequences.

> [!NOTE]
> Toyni was migrated from [jonas089's Github](https://github.com/jonas089/Toyni).
> Click [here](https://github.com/jonas089/Toyni) to see the past commit history.

## Status

**The Toyni STARK toolkit is in a solid state.** The core proof system
(domain construction, NTT-based FFT/IFFT, FRI low-degree testing, Merkle
commitment, Fiat-Shamir transcript, BabyBear field arithmetic) is
implemented end-to-end and exercised by the bundled Fibonacci AIR
(`src/fibonacci.rs` + `src/verifier.rs`). The unit-test suite passes and
the same primitives are used by the [zkvm](https://github.com/jonas089/zkvm)
project as its proving backend.

What is **not** in scope for Toyni: the AIR for any non-trivial program.
Toyni provides the building blocks; consumers (like zkvm) define their own
constraint systems on top.

## Learning

This implementation goes hand-in-hand with my article on STARKs that you
can read [here](https://github.com/jonas089/articles/blob/master/02-starks.md).
Toyni relies on prime fields and DEEP-ALI for soundness. I am looking
forward to exploring binary field STARKs later in my career.

![toyniii](art/toyniii.jpg)

*Meet the amazing artist behind this creation, [Kristiana Skrastina](https://www.linkedin.com/in/kristiana-skrastina/)*

## The Fibonacci example

Toyni ships with a small Fibonacci AIR that proves correct evaluation of
the recurrence over a power-of-two trace length:

```
| var |
|-----|
|  1  |
|  1  |
|  2  |
| ... |
| 13  |
| 21  |
```

The constraint reduces to:

```rust
fn fibonacci_constraint(t2: BabyBear, t1: BabyBear, t0: BabyBear) -> BabyBear {
    t2 - (t1 + t0)
}
```

Run it with:

```bash
cargo test test_fibonacci -- --nocapture
```

## CUDA NTT acceleration (`cuda` feature)

Toyni includes an optional CUDA backend for the NTT (forward + inverse).
It's gated behind the `cuda` feature flag and falls back to the CPU path
when the feature is off or no GPU is detected at runtime. The kernel and
its FFI live in [`cuda/ntt_kernel.cu`](cuda/ntt_kernel.cu) and
[`src/ntt.rs`](src/ntt.rs); `build.rs` invokes `nvcc` only when the
feature is enabled.

The path is tuned for repeated NTTs of the same size, the typical
proving workload. Per `n` it caches forward + inverse twiddles plus a
reusable device buffer in a global context, eliminating per-call
`cudaMalloc` / `cudaFree`, host-side twiddle precomputation, and per-stage
H2D copies. Mul reduction inside butterflies uses Barrett with
`mu = floor(2^64 / p)`. The build emits native code for sm_75/86/89 and,
when the toolkit supports it, sm_120 for Blackwell GPUs (RTX 50xx); a
forward-compat PTX target is always emitted as a fallback.

### NTT benchmark (RTX 5090, driver 595, CUDA 13)

End-to-end zkvm prover wall-time on the bundled `heavy-rust` example
(trace ≈ 2²⁰, LDE ≈ 2²³, ~491 NTTs of size 2²⁰ and 2²³ combined):

| Backend | Wall time |
|---------|-----------|
| CPU only | ~150 s |
| CUDA NTT (this branch) | ~11 s |

That's a ~13× speedup attributable almost entirely to the NTT path. The
remaining wall time is single-threaded CPU work in the consuming zkvm's
constraint and DEEP loops, which are not part of Toyni.

Build and test the CUDA path with:

```bash
cargo test --features cuda --release
```

## Theory: security properties

STARKs achieve their security through domain extension, low-degree testing,
and Merkle commitments:

1. **Domain Extension (Blowup).** The trace is extended to a domain `b`× larger
   than the trace length (here `b = 32`, sized to absorb the zero-knowledge
   masking below).
2. **Low-Degree Testing.** FRI folds the DEEP composition a *fixed* number of
   rounds down to a degree-bound layer, and the verifier reads that whole final
   layer and checks it is a constant (low-degree) codeword. **That final-layer
   check is what enforces the degree bound** — folding all the way to a single
   value and checking only that scalar enforces nothing and is forgeable.
3. **Merkle Commitments.** Each layer is committed via a Merkle tree; leaves are
   domain-separated, and the hiding (witness-carrying) trees are also salted.

The soundness error is roughly `ρ^q`, where `ρ` is the tested Reed–Solomon rate
and `q` the number of queries. The masked composition is tested at rate `1/8`
(degree bound `4·trace_len` over a `32·trace_len` domain), so 44 queries give
`~(1/8)^44 ≈ 2^-132`.

### Zero-knowledge

The trace polynomial is blinded as `T̂ = T + Z_H·R` with `R` uniformly random.
Because `Z_H` vanishes on the trace domain, `T̂ = T` there (constraints, and thus
soundness and completeness, are unchanged), but off it the LDE / query / OOD
openings are uniformly random. Together with per-leaf Merkle salting, the
verifier's view reveals nothing about the witness.

## Associated With

<div align="center">

| <a href="https://ciphercurve.com"><img src="https://ciphercurve.com/logo02.png" width="200" height="50" alt="Ciphercurve"></a> |
|:---:|
| [Ciphercurve](https://ciphercurve.com) |

</div>

---

<div align="center">
  <h3>2025 Ciphercurve</h3>
  <p><em>Building the future of privacy-preserving computation</em></p>
</div>
