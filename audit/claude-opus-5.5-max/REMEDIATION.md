# Remediation: toyni

Date: 2026-10-17. Goal: close the forgery path (TOYNI-01) with the smallest change to the pre-audit code. TOYNI-06 was fixed
afterwards with the same small change as the auth engine's AUTH-05. The proof format, the AIR and the `verify` signature are
unchanged.

## Patch

### TOYNI-01: the verifier rejects FRI rates above 1/8

`src/verifier.rs`
```diff
         let fri_degree_bound = (trace_len + MASK_DEGREE).next_power_of_two();
+        // Rate <= 1/8 (trace_len >= 64): at 8 and 16 any trace passes (TOYNI-01).
+        if lde_size < 8 * fri_degree_bound {
+            return false;
+        }
         let final_layer_size = lde_size / fri_degree_bound;
```

- With the current constants this rejects every `trace_len` below 64, which matches the rate the README states.
- `trace_len` still comes from the proof; the rate check is what bounds it.
- The check runs before `lde_size / final_layer_size`, so it also fixes the division by zero at `trace_len = 4` (TOYNI-04).

### TOYNI-06: queries are sampled over the whole LDE domain

`src/fibonacci.rs` and `src/verifier.rs`
```diff
-        let first_layer_half = fri_layers[0].len() / 2;          // prover; the verifier uses lde_size / 2
-        let query_indices = transcript.squeeze_indices(NUM_QUERIES, first_layer_half);
+        let query_indices = transcript.squeeze_indices(NUM_QUERIES, lde_size);
@@ prover
-            let deep_opening_pair = open_merkle_ext(&fri_trees[0], &fri_layers[0], qi + half0);
+            let deep_opening_pair = open_merkle_ext(&fri_trees[0], &fri_layers[0], qi ^ half0);
@@ verifier
-                    qi + half0,
+                    qi ^ half0,
-            let a0 = qp.deep_opening.value;
-            let b0 = qp.deep_opening_pair.value;
-            let x0_inv = shifted_elements[qi].inverse();
+            let lo0 = qi % half0;
+            let (a0, b0) = if qi == lo0 {
+                (qp.deep_opening.value, qp.deep_opening_pair.value)
+            } else {
+                (qp.deep_opening_pair.value, qp.deep_opening.value)
+            };
+            let x0_inv = shifted_elements[lo0].inverse();
-            let mut pos = qi;
+            let mut pos = lo0;
```

The DEEP check now covers the whole domain at rate 1/8. At `trace_len = 64` with 44 queries that is about 132 bits under the
conjectured bound and 66 under the proven one, which matches the README. `tests/adversarial.rs` and the PoC helpers (`common`,
`t05`) sample the same way.

## Status

| ID | Severity | Status |
|---|---|---|
| TOYNI-01 | Critical | Fixed |
| TOYNI-02 | High | Open in the library (`value` stays public). Closed at the auth engine's trust boundary (AUTH-02 patch). toyni's own verifier gets no forgery from it (REPORT). |
| TOYNI-03 | High | Open. Binding a statement needs public inputs, which means changing the `verify` API. The auth engine does not use this AIR. |
| TOYNI-04 | Medium | Fixed by the TOYNI-01 check |
| TOYNI-05 | Medium | Open (the Fibonacci example is not zero-knowledge; `t05` still recovers the witness, and the 140-coefficient mask is already short of the 144 equations from trace openings and OOD values, see REPORT Section 7) |
| TOYNI-06 | Medium | Fixed |
| TOYNI-07 – TOYNI-11 | Low / Info | Open |

## Re-audit

The PoCs were run against the pre-audit code and against the patched code (`cd poc && cargo test --release`, 12 tests).

| PoC | Pre-audit code | Patched |
|---|---|---|
| t01: random trace at `trace_len` 8 and 16 | verifies | rejected |
| t01: control at 64 | rejected | rejected |
| t02: `trace_len = 4` | panics | returns `false` |
| t03, t04, t05, t06 | demonstrate open findings | unchanged |

toyni's own suite passes (49 unit, 1 adversarial, 8 FRI).
