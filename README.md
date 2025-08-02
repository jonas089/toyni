# Toyni: A STARK Implementation in Progress
> [!WARNING]
> Toyni was migrated from [jonas089's Github](https://github.com/jonas089/Toyni)
> Click [here](https://github.com/jonas089/Toyni) to see the past commit history.

![toyniii](art/toyniii.jpg)

*Meet the amazing artist behind this creation, [Kristiana Skrastina](https://www.linkedin.com/in/kristiana-skrastina/)*

Toyni is an experimental zk-STARK proof system that is not meant for production use and has not been audited for correctness.
Currently Toyni is designed around a simple `Fibonacci` example program.

## Status Update
I'm currently working on a sound DEEP-ALI implementation in `prover.rs`. The degree of the deep polynomial D(x) blows up if C(x) does not vanish over
the fibonacci function (currently the only constraint) for the trace polynomial T(x). For soundness we will have to issue a challenge `alpha` and `z`
using fiat-shamir. This is tbd. The prover won’t be able to cheat because they must commit to a trace polynomial that results in a low-degree D(x) 
only if the constraints encoded in C(x) are satisfied.

The `Fibonacci` program defines a single-column trace table of shape:

```
| var |
|---|
| 1 |
| 1 |
| 2 | 
...  
| 13 |
| 21 | 
```

## Run the prover
```bash
cargo test test_fibonacci -- --nocapture
```

Output:

```bash
running 1 test
expected: 0, actual: 0
expected: 0, actual: 0
expected: 0, actual: 0
expected: 0, actual: 0
expected: 0, actual: 0
test prover::tests::test_fibonacci ... ok
```

Here the expected evaluation is ∑ci(x). For the fibonacci program we only have one constraint defined as:

```rust
  fn fibonacci_constraint(ti2: Fr, ti1: Fr, ti0: Fr) -> Fr {
      ti2 - (ti1 + ti0)
  }
```

> [!NOTE]
> We are skipping fiat shamir on purpose to get the core polynomial math right.
> Once the protocol itself is sound we will add fiat shamir where appropriate.
> For now we have fixed challenges and are trying to not leak the trace values.



# Theory: Security Properties
STARKs achieve their security through a combination of domain extension, low-degree testing, and Merkle commitments. Here's how it works:

The security of a STARK proof relies on three key mechanisms:

1. **Domain Extension (Blowup)**: The composition polynomial is evaluated over a domain that's `b` times larger than the original trace length.

2. **Low-Degree Testing**: The FRI protocol ensures that the polynomial being tested is close to a valid low-degree polynomial, with folding consistency checks at each layer.

3. **Merkle Commitments**: Each FRI layer is committed using a Merkle tree, ensuring the integrity of the folding process and enabling efficient verification.

The soundness error (probability of accepting an invalid proof) is bounded by:

```
Pr[undetected cheat] = (1/b)^q
```

where:
- `b` is the blowup factor (e.g., 8 in our example)
- `q` is the number of queries made by the verifier

This means that if a prover tries to cheat by modifying a fraction 1/b of the domain, the verifier will detect this with probability at least 1 - (1/b)^q. For example, with a blowup factor of 8 and 10 queries, the soundness error is at most (1/8)^10 ≈ 0.0000001.

# Associated With

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


