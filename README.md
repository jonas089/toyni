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

STARKs achieve their security through a combination of domain extension,
low-degree testing, and Merkle commitments. The three mechanisms are:

1. **Domain Extension (Blowup).** The composition polynomial is evaluated
   over a domain `b` times larger than the original trace length.
2. **Low-Degree Testing.** FRI ensures the polynomial being tested is
   close to a valid low-degree polynomial, with folding consistency checks
   at each layer.
3. **Merkle Commitments.** Each FRI layer is committed via a Merkle tree,
   ensuring the integrity of the folding process and enabling efficient
   verification.

The soundness error (probability of accepting an invalid proof) is
bounded by:

```
Pr[undetected cheat] = (1/b)^q
```

where `b` is the blowup factor and `q` is the number of queries. With
blowup 8 and 10 queries, the soundness error is at most
`(1/8)^10 ≈ 10⁻⁷`.

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
