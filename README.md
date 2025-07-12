# Toyni: A STARK Implementation in Progress
> [!WARNING]
> Toyni was migrated from [jonas089's Github](https://github.com/jonas089/Toyni)
> Click [here](https://github.com/jonas089/Toyni) to see the past commit history.

Welcome to Toyni! This is an implementation of a STARK (Scalable Transparent Argument of Knowledge) proving system in Rust.

![toyniii](art/toyniii.jpg)

Meet the amazing artist behind this creation, [Kristiana Skrastina](https://www.linkedin.com/in/kristiana-skrastina/)

> [!WARNING]  
> This is a research project and hasn't been audited. Use at your own risk.
> There are some essential and critical features still missing, especially the
> trace commitment with merkle proof verifications.
> This codebase is not complete by any means and an active work in progress.

## Work in Progress
Currently working on masked trace commitments for zero-knowledge spot checks. This is critical to the security of the protocol.
My understanding of this is:

1. interpolate the trace colums as polynomials
2. evaluate the trace column polynomials over the extended domain
3. interpolate the constraint polys from the original domain (tbd: random linear combinations) and sum them to get the composite constraint poly
4. check that C(X) = Z(X) * Q(X) (Q from the original domain) for X -> T(s) -> T + Z * r, where s is the spot check
from the extended domain

## 0. Background, STARK Verifier: Constraint vs FRI Layer Checks

## ✅ Constraint Check (Single Layer)
- For each randomly sampled point `x`:
  - Verifier checks:
    ```
    Q(x) * Z(x) == C(x)
    ```
  - Ensures that the execution trace satisfies all constraints
---

## ✅ FRI Layer Checks (Multiple Layers)
- Purpose: Prove that `Q(x)` is a **low-degree polynomial**
- Process:
  1. Start with evaluations of `Q(x)` over the domain (Layer 0)
  2. Recursively apply `fri_fold()` to reduce degree at each layer
  3. At each layer:
     - Verifier checks Merkle proofs for sampled values
     - Verifies that folding is consistent with previous layer
  4. Final layer should be constant or degree-1, checked directly

- ✅ Merkle proofs are **checked at each FRI layer**
- ✅ Folding correctness is verified at each layer

---

## 🔍 STARK Verifier Flow: Visual Diagram

```mermaid
sequenceDiagram
    participant Proof_Data
    participant Constraint_Check
    participant FRI_Protocol
    participant FRI_Merkle_Commitment
    participant FRI_Consistency
    participant FRI_Layers
    participant Accept_Proof

    Proof_Data-->Constraint_Check: Check Q(x)·Z(x) == C(x)
    Proof_Data-->FRI_Protocol: Begin FRI protocol
    FRI_Protocol-->FRI_Merkle_Commitment: Verify Merkle proofs
    FRI_Merkle_Commitment-->FRI_Consistency: Check folding inputs + β
    FRI_Consistency-->FRI_Layers: Fold with β0..βn
    FRI_Layers-->Accept_Proof: ✅ Accept Proof
```

>[!NOTE]
> FRI folding equation:
> f_{i+1}(x²) = (fᵢ(x) + fᵢ(–x) + βᵢ · (fᵢ(x) – fᵢ(–x))) / 2

## 🔐 STARK Verifier: Security Parameters for 128-bit Soundness

### FRI QUERIES
To achieve **128-bit soundness** in STARK proofs, the total probability that a cheating prover is accepted must be less than `2⁻¹²⁸`.

This involves carefully choosing parameters for:

- Constraint checks (`Q(x)` evaluations)
- FRI protocol (number of layers and queries per layer)

**Important**: The number of FRI layers depends on the program size!

> Example for different program sizes:
> - Trace size N=4: Extended domain = 32 → L = 4 layers → log₂(4) ≈ 2
> - Trace size N=8: Extended domain = 64 → L = 5 layers → log₂(5) ≈ 2.3  
> - Trace size N=16: Extended domain = 128 → L = 6 layers → log₂(6) ≈ 2.6
> - Trace size N=32: Extended domain = 256 → L = 7 layers → log₂(7) ≈ 2.8
> - Trace size N=1024: Extended domain = 8192 → L = 11 layers → log₂(11) ≈ 3.5
>
> **Formula**: FRI layers = log₂(8N) - 2, where N is trace size
> **Required queries**: m ≥ log₂(L) + 128 (e.g., 131-132 for most programs)

---

### CONSTRAINT CHECKS
> Example:
> - If `d / N = 1/4`, then `log₂(N/d) = 2`
> - So: `n = 128 / 2 = 64` spot checks

> - If `d / N = 1/8`, then `log₂(N/d) = 3`
> - So: `n = 128 / 3 ≈ 43`, but round up to be safe
---

### ✅ Practical Recommendation:
Use `n = 64–80` spot checks for strong 128-bit soundness across typical domain/degree ratios.

### ✅ Recommendations for 128-bit Security

| Component              | Suggested Value                    |
|-----------------------|------------------------------------|
| Constraint checks `n` | 64–80                              |
| FRI layers `L`        | log₂(8N) - 2 (where N = trace size) |
| FRI queries `m`       | ≥ log₂(L) + 128 (e.g., 131-132)     |
| Total soundness error | ε_total = ε_constraints + ε_fri ≤ 2⁻¹²⁸ |

**Note**: Our implementation uses a conservative approach by checking ALL points in each FRI layer rather than sampling, providing security well above the minimum requirements.

## 🔁 Summary

| Check Type         | Equation Checked              | Merkle Proofs | Multiple Layers? |
|--------------------|-------------------------------|----------------|-------------------|
| Constraint Check   | `Q(x) * Z(x) == C(x)`          | Optional       | ❌ No             |
| FRI Layer Check    | Folding consistency, low-degree| ✅ Yes          | ✅ Yes            |


## 1. Introduction

STARKs are a powerful cryptographic tool that enables proving the correct execution of a computation without revealing the underlying data. Think of it as a way to convince someone that you know the solution to a puzzle without actually showing them the solution. This property, known as zero-knowledge, is crucial for privacy-preserving applications in areas like financial transactions, voting systems, and private identity verification.

### 2. Why STARKs Matter

| Scalability | Transparency | Zero-Knowledge |
|-------------|--------------|----------------|
| • O(log² n) proof size | • No trusted setup | • Privacy |
| • Fast verify | • Public parameters | • Confidentiality |
| • Efficient | | • Data protection |
| | | • Secure sharing |

### 3. Real-World Applications

| Financial | Identity | Computing |
|-----------|----------|-----------|
| • Private payments | • Age verification | • Confidential computing |
| • Asset ownership | • Credential validation | • Private ML |
| | | • Secure MPC |

## 4. Technical Overview

At its heart, Toyni consists of three main components working together:

| Virtual Machine | Constraint System | STARK Prover |
|----------------|-------------------|--------------|
| •  | • Defines rules | • Generates proofs |
| • | • Validates states | • Uses FRI protocol |

### FRI Layer Scaling

The number of FRI layers in a STARK proof scales logarithmically with the program size:

```
FRI Layers = log₂(8N) - 2
```

Where `N` is the trace size (number of execution steps). This scaling ensures that:

- **Small programs** (N=4): 4 FRI layers
- **Medium programs** (N=32): 7 FRI layers  
- **Large programs** (N=1024): 11 FRI layers
- **Very large programs** (N=65536): 16 FRI layers

This logarithmic scaling is crucial for STARK's efficiency - proof size grows only logarithmically with computation size.

### 5. How It Works

| Program Execution | Execution Trace | Verification |
|------------------|-----------------|--------------|
| • Run program | • Record states | • Sample positions |
| • Track state | • Build constraints | • Check constraints |
| • Generate trace | | |

Here's a simple example that demonstrates how Toyni works. We'll create a program that proves a sequence of numbers increments by 1 each time:

```rust
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
```

This example demonstrates how Toyni can prove that a sequence of numbers follows a specific pattern (incrementing by 1) without revealing the actual numbers. The proof can be verified by anyone, but the actual values remain private.

### 6. Security Properties

STARKs achieve their security through a combination of domain extension and low-degree testing. Here's how it works:

| Domain Extension | Low-Degree Testing | Soundness Guarantees |
|-----------------|-------------------|---------------------|
| • Extend domain | • FRI protocol | • Soundness error: (1/b)^q |
| • Blowup factor | • Polynomial degree | • Query complexity |

The security of a STARK proof relies on two key mechanisms:

1. **Domain Extension (Blowup)**: The composition polynomial is evaluated over a domain that's `b` times larger than the original trace length, where `b` is the blowup factor (8 in our implementation).

2. **Low-Degree Testing**: The FRI protocol ensures that the polynomial being tested is close to a valid low-degree polynomial. The number of FRI layers scales logarithmically with program size: `L = log₂(8N) - 2` where N is the trace size.

The soundness error (probability of accepting an invalid proof) is bounded by:

```
Pr[undetected cheat] = (1/b)^q
```

where:
- `b` is the blowup factor (8 in our implementation)
- `q` is the number of queries made by the verifier

This means that if a prover tries to cheat by modifying a fraction 1/b of the domain, the verifier will detect this with probability at least 1 - (1/b)^q. 

**Example**: With a blowup factor of 8 and 80 constraint queries, the constraint soundness error is at most (1/8)^80 ≈ 2^(-240), providing far more than the required 128-bit security.

## 7. Project Structure

The codebase is organized into logical components:

| Math | VM | Library |
|------|----|---------|
| • Polynomial | • Constraints | • Entry point |
| • Domain | • Trace | • Public API |
| • FRI | • Execution | • Documentation |
| • STARK | | |


### 8. Current Features

| Constraint System | FRI Protocol | Mathematical Operations |
|------------------|--------------|------------------------|
| • Transition constraints | • Low-degree testing | • Polynomial arithmetic |
| • Boundary constraints | • Interactive verification | • Field operations |
| • Quotient verification | • FRI folding layers | • Domain operations |
| • Merkle commitments | • Folding consistency checks | • Secure commitments |
| • Trace Privacy | • Fiat Shamir verifier challenges | • Conservative FRI verification |

**Security Note**: Our FRI implementation uses a conservative approach by verifying ALL points in each layer rather than sampling, providing security well above the theoretical minimum requirements.

### 9. Missing Components

| Zero-Knowledge | Performance | Optimization |
|----------------|-------------|--------------|
| • Enhanced state protection | • GPU Acceleration | • IFFT for interpolation |
| • Deterministic hashing | • Non-interactive proofs | • No dependency on arkworks |

While we have a working STARK implementation with quotient polynomial verification, FRI folding, and basic zero-knowledge properties, there are still some components to implement:

1. **Performance Optimizations**: Need to implement parallel processing and batch verification for better scalability.
2. **Circuit-Specific Features**: Add support for specialized circuits and optimizations.

### 10. Roadmap

#### Completed Features ✅
- Basic STARK implementation with constraint checks
- FRI protocol with folding layers
- Merkle commitments for FRI layers
- Folding consistency verification
- Interactive verification protocol
- Fiat-Shamir transform implementation
- Zero-knowledge polynomial masking

#### In Progress 🚧
- Performance optimizations
- Circuit-specific optimizations

#### Future Work 📅
- Enhanced zero-knowledge properties
- Parallel processing support
- Batch verification
- Circuit-specific optimizations
- Documentation improvements

### 11. Security Properties

STARKs achieve their security through a combination of domain extension, low-degree testing, and Merkle commitments. Here's how it works:

| Domain Extension | Low-Degree Testing | Merkle Commitments |
|-----------------|-------------------|-------------------|
| • Extend domain | • FRI protocol | • Tree structure |
| • Blowup factor | • Polynomial degree | • Proof generation |
| • Soundness | • Folding checks | • Commitment verification |

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

## 12. Contributing

We welcome contributions to Toyni! Our current focus is on enhancing zero-knowledge properties and improving the overall system. We're particularly interested in:

1. Enhancing zero-knowledge properties and state protection
2. Adding comprehensive test coverage and security audits
3. Improving documentation and adding more examples
4. Optimizing performance and reducing proof sizes

# 13. Associated With

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


