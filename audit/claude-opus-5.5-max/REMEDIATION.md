# Remediation: toyni

Date: 2026-10-17. Goal: close the forgery path (TOYNI-01) with the smallest change to the pre-audit code. The proof format, the
AIR and the `verify` signature are unchanged.

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

## Status

| ID | Severity | Status |
|---|---|---|
| TOYNI-01 | Critical | Fixed |
| TOYNI-02 | High | Open in the library (`value` stays public). Closed at the auth engine's trust boundary (AUTH-02 patch). toyni's own verifier gets no forgery from it (REPORT). |
| TOYNI-03 | High | Open. Binding a statement needs public inputs, which means changing the `verify` API. The auth engine does not use this AIR. |
| TOYNI-04 | Medium | Fixed by the TOYNI-01 check |
| TOYNI-05, TOYNI-06 | Medium | Open (TOYNI-06: ≈88 bits conjectured, ≈44 bits proven at `trace_len = 64`) |
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
