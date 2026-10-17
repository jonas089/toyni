# Security Audit: toyni STARK toolkit

| | |
|---|---|
| Auditor | Claude Opus 5.5 (AI model, max reasoning effort) |
| Date | 2026-10-17 |
| Repository | `toyni`, HEAD `af8683c` plus the working tree (modified: `src/ext.rs` norm-based inverse, `src/ntt/cpu.rs` in-place NTT) |
| PoCs | `audit/claude-opus-5.5-max/poc` (standalone crate; `cargo test --release`) |
| Consumer | `toyni-zkvm` auth engine; see `toyni-zkvm/audit/claude-opus-5.5-max/REPORT.md` |

> This is an AI-generated audit. It is not a substitute for a review by qualified human cryptographers. Every finding marked
> *Confirmed* has an executable proof of concept that passes against the audited code.

## 1. Summary

toyni provides BabyBear and its quartic extension, NTT/FFT, Merkle trees, a Fiat-Shamir transcript, FRI folding and masking, plus a
Fibonacci STARK (prover and verifier) as the reference AIR.

- **Forgery.** The reference verifier accepts proofs over traces that violate the AIR when the prover picks a small trace length
  (TOYNI-01).
- **Field arithmetic.** Arithmetic and equality are wrong for non-canonical field elements. That is the root cause of a complete forgery
  in the auth engine (TOYNI-02 → AUTH-02).
- **Example AIR and README claims.** The Fibonacci AIR binds no statement (TOYNI-03), and its zero-knowledge claim is false (TOYNI-05).

| ID | Title | Severity | Status |
|---|---|---|---|
| TOYNI-01 | Verifier trusts prover-chosen `trace_len`; at 8 and 16 the low-degree test is vacuous | **Critical** | Confirmed (PoC) |
| TOYNI-02 | `BabyBear` arithmetic and equality are wrong for non-canonical values (public `u32`) | **High** (Critical in consumers) | Confirmed (PoC; exploited in AUTH-02) |
| TOYNI-03 | Fibonacci AIR has no boundary constraints or public inputs | **High** | Confirmed (PoC) |
| TOYNI-04 | Verifier panics (division by zero) on `trace_len = 4` | Medium | Confirmed (PoC) |
| TOYNI-05 | Zero-knowledge claim false: the witness is recovered from one proof | Medium | Confirmed (PoC) |
| TOYNI-06 | Stated 132-bit soundness is ~88 bits conjectured and ~44 bits proven | Medium | Analysis |
| TOYNI-07 | Fiat-Shamir initialization binds no domain, statement or parameters | Low | Analysis |
| TOYNI-08 | Merkle root does not bind the leaf count (CVE-2012-2459 pattern); unchecked sibling length | Low | Confirmed (PoC) |
| TOYNI-09 | CUDA NTT host code ignores API errors and over-allocates for large `n` | Low | Analysis |
| TOYNI-10 | `squeeze_indices` loops forever when `count > max` | Info | Analysis |
| TOYNI-11 | Incorrect or dead helpers (`divide_by_linear`, truncating `fft`) | Info | Analysis |

Section 6 lists the properties that were checked and hold.

## 2. Scope

All of `src/` (babybear, ext, merkle, transcript, fibonacci, verifier, math/{domain, fri, mask, polynomial}, ntt/{cpu, cuda},
program/trace), `cuda/ntt_kernel.cu`, `build.rs`, `tests/` and `README.md`, about 3,900 lines.
The two working-tree modifications were audited as modified, and both are covered by tests against reference implementations.

## 3. Threat model and methodology

The adversary is a malicious prover that controls every byte of a `StarkProof`, including `trace_len` and `lde_size`, and anyone who
submits proofs to a verifier. The methodology and the published vulnerability classes are the same as in the auth engine report:
- the Least Authority Plonky3 audit (2024): Issue A on parameters trusted from the proof, Issues B and C on Fiat-Shamir, Issues D–F on
  panics, and Suggestion 4 on soundness accounting;
- Frozen Heart and weak Fiat-Shamir (Dao et al., S&P 2023);
- under-constrained circuits (Chaliasos et al., USENIX Security 2024; Tornado Cash's 2019 MiMC circuit bug);
- STARK zero-knowledge (Haböck and Al Kindi, ePrint 2024/1037);
- CVE-2012-2459 (Bitcoin Merkle duplication).

Severity levels are as defined in the auth engine report.

## 4. Findings

### TOYNI-01: The verifier trusts the prover-chosen `trace_len`; at 8 and 16 the low-degree test is vacuous

**Severity: Critical.** Location: `src/verifier.rs:25-36, 92-99, 126-127, 163-194`; `src/fibonacci.rs:246-247`.

**Description.** The verifier reads `trace_len` from the proof and only requires a power of two in `[4, 2^20]`.

| `trace_len` | `lde = 32n` | `D = next_pow2(n + 140)` | Effect |
|---|---|---|---|
| 8 | 256 | 256 | Rate 1: every word is a codeword; the final layer has size 1 and is trivially "constant" |
| 16 | 512 | 256 | Rate 1/2, but DEEP consistency is checked only at query indices in `[0, lde/2)` (lines 126-127, 179-192), where `D` equals the number of points |
| 4 | 128 | 256 | Division by zero (TOYNI-04) |

At 8 and 16, any committed trace and quotient pass: claim arbitrary OOD values, solve `q(z)` from the constraint identity, and interpolate
the DEEP layer from the first half.

**Exploit** (`poc/tests/t01_short_trace_forgery.rs`). A uniformly random, non-Fibonacci trace verifies at `trace_len = 8` and `16`. The
same malicious prover is rejected at 64 (control test).

**Impact.** The verifier accepts proofs whose committed trace violates the AIR. With the Fibonacci AIR the statement is vacuous anyway
(TOYNI-03). But consumers inherit the pattern: the auth engine's first-half DEEP check is exactly AUTH-01, a complete forgery.

**Recommendation.**
- The verifier must take `trace_len` from the statement, not the proof (Least Authority Plonky3 Issue A).
- Enforce a maximum rate, e.g. `D <= lde/4`.
- Check DEEP consistency at both positions of every FRI pair, or sample queries over the full domain.

### TOYNI-02: `BabyBear` arithmetic and equality are wrong for non-canonical values

**Severity: High** (library root cause; Critical in the auth engine, AUTH-02/03). Location: `src/babybear.rs:14-17, 20-24, 62-66, 90-92,
136-160, 202-210`; `src/ext.rs:22-24, 92-98`.

**Description.** `pub struct BabyBear { pub value: u32 }` documents "always canonical" but enforces nothing.
- `+`: `self.value + rhs.value` wraps the `u32` (release) or panics (debug) when the sum is at least `2^32`.
- `-`: `self.value + PRIME - rhs.value` does the same when `rhs.value >= p` or `self.value` is large.
- `neg`: underflows for `value > p`.
- Each wrap shifts the result by `Δ = 2^32 mod p`.
- `PartialEq` and `Hash` compare raw values, so `BabyBear{value: p} != 0`.
- `from_bytes` silently reduces, so encodings are malleable. `Ext` inherits all of this.

Multiplication reduces correctly, which hides the problem in most code paths.

**Exploit.**
- `poc/tests/t04_noncanonical_arithmetic.rs`: equality, malleable decoding, a wrong sum, and the same OOD value producing different field
  elements in a DEEP numerator and a constraint expression.
- The auth engine PoC `f02_noncanonical_ood_forgery.rs` turns this into a complete forgery.
- In toyni's own Fibonacci verifier the attacker gets only 12 binary shifts (3 OOD values plus `q_z`, coordinates 1–3) against a 124-bit
  target, and the statement is vacuous, so no direct forgery there.

**Recommendation.**
- Make the field private, with `new` (reducing) and `try_from_canonical` (rejecting) constructors.
- Decode with an error for `value >= p`.
- Implement add, sub and neg with `u64` intermediates, or `debug_assert` canonical inputs.
- Derive equality from the canonical value.

### TOYNI-03: The Fibonacci AIR has no boundary constraints or public inputs

**Severity: High.** Location: `src/fibonacci.rs:139-151, 335-345`; `src/verifier.rs:54-74`; `README.md` ("The Fibonacci example").

**Description.**
- The README presents the AIR as illustrating "transition and boundary constraints", but the code has none of the latter.
  `boundary_constraint_1/2` are exemption factors, `(x - g^(n-1))(x - g^(n-2))`, that switch off the two wrap-around transitions.
- No public input (seeds, length or result) is bound or checked, and the transcript absorbs nothing but commitments.
- A proof therefore attests only that *some* sequence satisfies `t[i+2] = t[i+1] + t[i]` on rows 0..n−3. That is always true.

**Exploit** (`poc/tests/t03_unbound_statement.rs`). The all-zero sequence and a sequence with seeds `(1000000, 7)` both verify.

**Impact.** The proof cannot be tied to `fib(n)`, which is the classic under-constrained-circuit failure.

**Recommendation.** Add public inputs (`t[0]`, `t[1]`, the claimed `t[n-1]`), absorb them first in the transcript, and enforce them as
boundary quotients, e.g. `(T(x) - 1)/(x - 1)`, checked at `z` like the transition quotient.

### TOYNI-04: The verifier panics (division by zero) on `trace_len = 4`

**Severity: Medium** (verifier DoS). Location: `src/verifier.rs:29, 92-94`.

**Description.** `D = 256 > lde = 128`, so `final_layer_size = 0` and `lde_size / final_layer_size` panics. The verifier's own doc comment
promises that "every failure path returns false". Only all-zero OOD values and 32-byte commitments are needed.

**Exploit.** `poc/tests/t02_trace_len_4_panics.rs`.

**Recommendation.** Reject `D > lde/2` (better, `D > lde/4`) before dividing, and fuzz `verify` with malformed proofs (Least Authority
Plonky3 Issues E and F).

### TOYNI-05: The zero-knowledge claim is false; the witness is recovered from one proof

**Severity: Medium** (privacy; the README claims "the verifier's view reveals nothing about the witness"). Location:
`src/fibonacci.rs:18-20, 127-133, 197-224`.

**Description.**
- `MASK_DEGREE = 3·44 + 8 = 140` is budgeted for "3 openings per query + 3 OOD points". That makes `T^ = T + Z_H R` a polynomial with 204
  unknown coefficients.
- A proof reveals more linear functionals of `T^` than that:
  - 132 opened trace values,
  - 3 OOD values (12 base equations),
  - 44 DEEP pair values at `qi + lde/2`, 4 base equations each. These are not counted, and they are linear in `T^` at three points because
    this quotient is a fixed linear function of `T^`.
  - FRI layer values add more.

**Exploit** (`poc/tests/t05_zk_witness_recovery.rs`). It builds the 320 × 204 system from public proof data only, solves it by Gaussian
elimination, and recovers the secret trace exactly.

**Recommendation.** Budget masking for every revealed functional, including the DEEP pair and FRI layers, or randomize the DEEP/FRI input.
Also mask the quotient (see Haböck and Al Kindi, ePrint 2024/1037). Until then, remove the zero-knowledge claim.

### TOYNI-06: The stated 132-bit soundness is ~88 bits conjectured and ~44 bits proven

**Severity: Medium.** Location: `src/fibonacci.rs:11-15`; `README.md` ("soundness error is roughly ρ^q … ≈ 2^-132").

**Description.**
- With DEEP checked only on the first half (TOYNI-01), the effective rate at `trace_len = 64` is 1/4, not 1/8.
- 44 queries then give about 88 bits under the conjectured bound and about 44 bits under the proven (Johnson-regime) bound.
- There is no grinding.
- For prover-chosen `trace_len` of 8 or 16 the security is 0 (TOYNI-01).

**Recommendation.** Fix TOYNI-01, then compute and document conjectured and proven soundness from all parameters (field size, rate,
queries, grinding, DEEP/ALI terms; Least Authority Plonky3 Suggestion 4; ethSTARK ePrint 2021/582).

### TOYNI-07: Fiat-Shamir initialization binds no domain, statement or parameters

**Severity: Low.** Location: `src/transcript.rs:12-56`; `src/verifier.rs:54-66`.

**Description.**
- The transcript starts from a fixed label and absorbs only commitments and OOD values. It binds no protocol or AIR identifier, no public
  inputs (none exist; TOYNI-03), and no `trace_len`, blowup, query count or mask degree.
- Messages are concatenated without type or length framing, and `absorb_commitment` accepts any length. This verifier checks 32 bytes;
  consumers must remember to.
- This is Least Authority's Plonky3 Issue B, and the root of the weak Fiat-Shamir/Frozen Heart class.

**Recommendation.** Provide a transcript API that requires a domain separator, a parameter digest and the statement before the first
commitment, and frames every message.

### TOYNI-08: The Merkle root does not bind the leaf count, and sibling length is unchecked

**Severity: Low.** Location: `src/merkle.rs:41-54, 88-123`.

**Description.**
- An odd level duplicates its last node, so `[a, b, c]` and `[a, b, c, c]` share a root (the CVE-2012-2459 pattern), and an opening at
  index 3 verifies against the 3-leaf root when the verifier is told 4 leaves.
- Binding therefore relies on every caller fixing `num_leaves`, which the current STARKs do (powers of two).
- `verify_merkle_proof` also accepts siblings of any length (`verify_multi_proof` checks for 32).

**Exploit.** `poc/tests/t06_merkle_duplicate_last_leaf.rs`.

**Recommendation.** Hash the leaf count into the root or reject odd levels, and require 32-byte siblings.

### TOYNI-09: The CUDA NTT host code ignores API errors and over-allocates for large `n`

**Severity: Low** (completeness and availability). Location: `cuda/ntt_kernel.cu:67, 112-121`; `src/ntt/cuda.rs`.

**Description.**
- `cudaMalloc`, `cudaMemcpy` and stream calls are unchecked. Only `cudaGetLastError` is returned at the end, so failures produce garbage
  commitments or late errors.
- Each of 3 streams allocates `CHUNK_POLYS (1024) × n × 4` bytes. That is 48 GB at `n = 2^22`, which fails on any consumer GPU.
- `__brev(i) >> (32 - log_n)` is undefined for `n = 1`.
- The README describes Barrett reduction and per-`n` twiddle and buffer caching, which the kernel does not implement.

**Recommendation.** Check every CUDA call, size chunks from `n` and the available memory, special-case `n = 1`, and update the README.

### TOYNI-10: `squeeze_indices` loops forever when `count > max`

**Severity: Info.** Location: `src/transcript.rs:58-72`.

The loop draws distinct indices until it has `count` of them. This is unreachable with the shipped parameters, but a consumer with a small
domain would hang. Assert `count <= max`.

### TOYNI-11: Incorrect or dead helpers

**Severity: Info.**
- `Polynomial::divide_by_linear` (`src/math/polynomial.rs:70-85`) breaks on its first iteration and returns a zero quotient and remainder
  for every input. It is unused.
- `BabyBearDomain::fft` (`src/math/domain.rs:116-118`) silently truncates coefficient vectors longer than the domain.
- `ExecutionTrace::interpolate_column` is O(n²) and unused.

## 5. Proof-of-concept index

| Finding | Test |
|---|---|
| TOYNI-01 | `tests/t01_short_trace_forgery.rs` (`trace_len` 8 and 16 accepted; 64 rejected as control) |
| TOYNI-02 | `tests/t04_noncanonical_arithmetic.rs`; full forgery in `toyni-zkvm/audit/.../f02_noncanonical_ood_forgery.rs` |
| TOYNI-03 | `tests/t03_unbound_statement.rs` |
| TOYNI-04 | `tests/t02_trace_len_4_panics.rs` |
| TOYNI-05 | `tests/t05_zk_witness_recovery.rs` |
| TOYNI-08 | `tests/t06_merkle_duplicate_last_leaf.rs` |

```
cd audit/claude-opus-5.5-max/poc && cargo test --release     # 12 tests
cargo test                                                    # debug: wrap tests are skipped
```

## 6. Properties checked that hold

- **Field and domain.**
  - `440564289 = 31^15` is a primitive `2^27`-th root of unity.
  - The coset shift 7 lies outside every `2^k` subgroup.
  - `X^4 - 11` is irreducible over BabyBear (11 is a quadratic non-residue and `p ≡ 1 mod 4`).
  - NTT roots and domain roots coincide.
- **Ext inverse.** The norm-based inverse agrees with Fermat inversion (`norm_inverse_matches_fermat`).
- **NTT.** The in-place NTT agrees bit-for-bit with the reference `ntt`/`intt` (`babybear_matches_reference_ntt`).
- **Merkle.**
  - Openings are bound to the verifier's index, with path depth pinned.
  - Leaves and nodes are domain-separated.
  - The multi-proof verifier rejects duplicates, slack and moved leaves.
- **Transcript order.** The final FRI layer is committed through its root before queries (Least Authority Plonky3 Issue C does not apply).
  DEEP uses independent per-term coefficients drawn after the OOD claims, which fixes the earlier forgery in `tests/adversarial.rs`.
- **Verifier hardening.** `derive_z_verifier` is bounded, commitment lengths are checked, and the parameter shape is validated (apart from
  TOYNI-01 and TOYNI-04).

## 7. Re-audit after remediation (2026-10-17)

**Scope.** Critical and High only: any way for a verifier to accept a false statement. The code is the audited code plus the
TOYNI-01 patch in `REMEDIATION.md`.

**Result: no new Critical or High findings.**

| ID | Severity | Status |
|---|---|---|
| TOYNI-01 | Critical | Fixed (t01 rejected at `trace_len` 8 and 16) |
| TOYNI-02 | High | Open in the library; closed at the auth engine's boundary (AUTH-02). toyni's own verifier still accepts values >= p. That has no impact while TOYNI-03 leaves the statement vacuous, so fix it before or together with TOYNI-03. |
| TOYNI-03 | High | Open. This is not a verifier soundness break: the proof is valid for a statement that binds nothing. The auth engine does not use this AIR. |

**Method.**
- The Fibonacci verifier was reviewed line by line: parameter checks, transcript, OOD check, DEEP terms, FRI and openings.
- The primitives shared with the auth engine (field, extension, Merkle, transcript, domain) were reviewed twice, independently;
  see the auth engine report, Section 7.

## Appendix A: Soundness estimates (`trace_len = 64`, 44 queries, no grinding)

| | ρ | Conjectured | Proven (Johnson) |
|---|---|---|---|
| README claim | 1/8 | 132 | 66 |
| As implemented (DEEP on first half) | 1/4 | ≈ 88 | ≈ 44 |
| `trace_len` 8 or 16 | 1 | 0 | 0 |

## References

1. Least Authority, *Plonky3 Security Audit Report* (Polygon), 7 Nov 2024.
2. Q. Dao, J. Miller, O. Wright, P. Grubbs, *Weak Fiat-Shamir Attacks on Modern Proof Systems*, IEEE S&P 2023, ePrint 2023/691.
3. Trail of Bits, *Frozen Heart* disclosures, 2022.
4. S. Chaliasos et al., *SoK: What Don't We Know? Understanding Security Vulnerabilities in SNARKs*, USENIX Security 2024, arXiv 2402.15293.
5. U. Haböck, Al Kindi, *A note on adding zero-knowledge to STARKs*, ePrint 2024/1037.
6. StarkWare, *ethSTARK Documentation*, ePrint 2021/582.
7. E. Ben-Sasson et al., *DEEP-FRI*, ePrint 2019/336; *Proximity Gaps for Reed–Solomon Codes*, ePrint 2020/654.
8. CVE-2012-2459, Bitcoin Merkle tree duplicate-leaf malleability.
9. Tornado Cash, October 2019 disclosure of the under-constrained MiMC circuit (forged deposit proofs).
