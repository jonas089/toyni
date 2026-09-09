use crate::babybear::BabyBear;
use crate::ext::Ext;
use crate::math::domain::BabyBearDomain;
use crate::merkle::{verify_merkle_proof, MerkleTree};
use crate::fibonacci::{
    deep_value, eval_boundary_1_ext, eval_boundary_2_ext, eval_fibonacci_constraint_ext,
    MerkleOpening, MerkleOpeningExt, StarkProof, BLOWUP, COSET_SHIFT, MASK_DEGREE,
    NUM_DEEP_TERMS, NUM_QUERIES,
};
use crate::transcript::FiatShamirTranscript;

/// Largest trace length the verifier will consider. Bounds the work a malformed
/// proof can make it do, and keeps the LDE inside BabyBear's 2-adicity (2^27).
pub const MAX_TRACE_LEN: usize = 1 << 20;

pub struct StarkVerifier;

impl StarkVerifier {
    /// Check `proof`. A proof is untrusted input, so every failure path returns
    /// `false`: an earlier version reached `assert!`s inside
    /// `BabyBearDomain::new` and `get_root_of_unity`, indexed
    /// `fri_final_layer` unchecked, and could spin in `derive_z_verifier` —
    /// all reachable from a malformed proof.
    pub fn verify(&self, proof: &StarkProof) -> bool {
        let trace_len = proof.trace_len;
        let lde_size = proof.lde_size;

        // Parameters first: everything below indexes and allocates from these.
        if !trace_len.is_power_of_two() || trace_len < 4 || trace_len > MAX_TRACE_LEN {
            return false;
        }
        if lde_size != trace_len * BLOWUP
            || !lde_size.is_power_of_two()
            || lde_size.trailing_zeros() > 27
        {
            return false;
        }
        if proof.trace_commitment.len() != 32
            || proof.quotient_commitment.len() != 32
            || proof.fri_commitments.iter().any(|c| c.len() != 32)
        {
            return false;
        }

        let domain = BabyBearDomain::new(trace_len);
        let extended_domain = BabyBearDomain::new(lde_size);
        let shift = BabyBear::new(COSET_SHIFT);
        let shifted_domain = extended_domain.get_coset(shift);
        let g = domain.group_gen();

        let vanishing = domain.vanishing_poly_coeffs();

        // ── 1. Replay Fiat-Shamir transcript ───────────────────────────
        let mut transcript = FiatShamirTranscript::new();
        transcript.absorb_commitment(&proof.trace_commitment);
        transcript.absorb_commitment(&proof.quotient_commitment);

        let z = match derive_z_verifier(&mut transcript) {
            Some(z) => z,
            None => return false,
        };

        transcript.absorb_ext(proof.t_z);
        transcript.absorb_ext(proof.t_gz);
        transcript.absorb_ext(proof.t_ggz);
        transcript.absorb_ext(proof.q_z);

        // ── 2. OOD constraint check: C(z) = Q(z) · Z(z) ──────────────
        let c_z = eval_fibonacci_constraint_ext(proof.t_ggz, proof.t_gz, proof.t_z)
            * eval_boundary_1_ext(z, g, trace_len)
            * eval_boundary_2_ext(z, g, trace_len);
        if c_z != proof.q_z * z.eval_base_at_ext(&vanishing) {
            return false;
        }

        // ── 3. DEEP coefficients ───────────────────────────────────────
        // Drawn after the OOD values are absorbed, so the prover commits to its
        // claims before learning how they will be weighted.
        let mut deep_coeffs = [Ext::zero(); NUM_DEEP_TERMS];
        for c in deep_coeffs.iter_mut() {
            *c = transcript.squeeze_ext_challenge();
        }

        // ── 4. Replay FRI commitments & derive betas ───────────────────
        if proof.fri_commitments.is_empty() {
            return false;
        }

        // Fold down to the degree-bound layer (size lde/D_BOUND) and require it
        // to be constant. This enforces the degree bound; fold-consistency alone
        // does not. D_BOUND = next_pow2(trace_len + MASK_DEGREE) covers masking.
        let fri_degree_bound = (trace_len + MASK_DEGREE).next_power_of_two();
        let final_layer_size = lde_size / fri_degree_bound;
        let expected_folds = (lde_size / final_layer_size).trailing_zeros() as usize;
        if proof.fri_commitments.len() != expected_folds + 1 {
            return false;
        }
        if proof.fri_final_layer.len() != final_layer_size {
            return false;
        }
        // Final layer is constant (degree 0).
        if !proof
            .fri_final_layer
            .iter()
            .all(|v| *v == proof.fri_final_layer[0])
        {
            return false;
        }
        // Final layer matches its commitment (binds all positions to the transcript).
        if merkle_root_of_ext(&proof.fri_final_layer) != *proof.fri_commitments.last().unwrap() {
            return false;
        }

        transcript.absorb_commitment(&proof.fri_commitments[0]);

        let num_fri_folds = proof.fri_commitments.len() - 1;
        let mut fri_betas: Vec<Ext> = Vec::with_capacity(num_fri_folds);

        for i in 1..proof.fri_commitments.len() {
            let beta = transcript.squeeze_ext_challenge();
            fri_betas.push(beta);
            transcript.absorb_commitment(&proof.fri_commitments[i]);
        }

        // ── 4. Derive query indices ────────────────────────────────────
        let first_layer_half = lde_size / 2;
        let query_indices = transcript.squeeze_indices(NUM_QUERIES, first_layer_half);

        if proof.query_proofs.len() != NUM_QUERIES {
            return false;
        }

        // ── 5. Shifted domain elements (for x-coordinate lookups) ──────
        let shifted_elements = shifted_domain.elements();
        let half_inv = BabyBear::new(2).inverse();

        // ── 6. Verify each query ───────────────────────────────────────
        for (qi_idx, qp) in proof.query_proofs.iter().enumerate() {
            let qi = query_indices[qi_idx];
            if qp.index != qi {
                return false;
            }
            // The query must carry exactly the intermediate-layer openings the
            // fixed fold schedule produces (layers 1..final, exclusive).
            if qp.fri_openings.len() != expected_folds - 1 {
                return false;
            }

            // 6a. Trace and quotient openings. Each is verified against the
            //     index the verifier derived, not one carried in the proof.
            let idx_g = (qi + BLOWUP) % lde_size;
            let idx_gg = (qi + 2 * BLOWUP) % lde_size;
            if !verify_opening(&qp.trace_opening, qi, lde_size, &proof.trace_commitment)
                || !verify_opening(&qp.trace_opening_g, idx_g, lde_size, &proof.trace_commitment)
                || !verify_opening(&qp.trace_opening_gg, idx_gg, lde_size, &proof.trace_commitment)
            {
                return false;
            }
            if !verify_opening(&qp.quotient_opening, qi, lde_size, &proof.quotient_commitment) {
                return false;
            }

            // 6b. DEEP layer (FRI layer 0): the position and its fold pair.
            let half0 = lde_size / 2;
            if !verify_opening_ext(&qp.deep_opening, qi, lde_size, &proof.fri_commitments[0])
                || !verify_opening_ext(
                    &qp.deep_opening_pair,
                    qi + half0,
                    lde_size,
                    &proof.fri_commitments[0],
                )
            {
                return false;
            }

            // 6c. DEEP consistency: the committed columns really do compose
            //     into the committed DEEP layer, *term by term*.
            let x_i = shifted_elements[qi];
            let expected_deep = deep_value(
                Ext::from(x_i),
                z,
                Ext::from(qp.trace_opening.value),
                Ext::from(qp.trace_opening_g.value),
                Ext::from(qp.trace_opening_gg.value),
                Ext::from(qp.quotient_opening.value),
                proof.t_z,
                proof.t_gz,
                proof.t_ggz,
                proof.q_z,
                &deep_coeffs,
            );
            if qp.deep_opening.value != expected_deep {
                return false;
            }

            // 6e. First FRI fold: layer 0 → layer 1
            let a0 = qp.deep_opening.value;
            let b0 = qp.deep_opening_pair.value;
            let x0_inv = shifted_elements[qi].inverse();

            let mut prev_folded = {
                let avg = (a0 + b0).mul_base(half_inv);
                let diff = (a0 - b0).mul_base(half_inv);
                avg + diff * fri_betas[0] * Ext::from(x0_inv)
            };

            // 6f. Intermediate FRI layers
            let mut pos = qi;

            for layer in 0..qp.fri_openings.len() {
                let fold_k = layer + 1;
                let layer_size = lde_size >> fold_k;
                let half = layer_size / 2;

                if half == 0 || pos >= layer_size {
                    return false;
                }
                let lo = pos % half;
                let in_first_half = pos == lo;

                let (ref op, ref op_pair) = qp.fri_openings[layer];

                // Merkle proofs, bound to the derived positions.
                if !verify_opening_ext(op, lo, layer_size, &proof.fri_commitments[fold_k])
                    || !verify_opening_ext(
                        op_pair,
                        lo + half,
                        layer_size,
                        &proof.fri_commitments[fold_k],
                    )
                {
                    return false;
                }

                // prev_folded should match the value at position `pos`
                if in_first_half {
                    if op.value != prev_folded {
                        return false;
                    }
                } else if op_pair.value != prev_folded {
                    return false;
                }

                // x-coordinate: xs_k[lo] = shifted_elements[lo]^{2^fold_k}
                let x_inv = shifted_elements[lo].pow(1u64 << fold_k).inverse();

                let a_l = op.value;
                let b_l = op_pair.value;
                let avg = (a_l + b_l).mul_base(half_inv);
                let diff = (a_l - b_l).mul_base(half_inv);
                prev_folded = avg + diff * fri_betas[fold_k] * Ext::from(x_inv);

                pos = lo;
            }

            // 6g. The folded query position must land on the final layer.
            if pos >= proof.fri_final_layer.len() || proof.fri_final_layer[pos] != prev_folded {
                return false;
            }
        }

        true
    }
}

fn verify_opening(
    opening: &MerkleOpening,
    index: usize,
    num_leaves: usize,
    root: &[u8],
) -> bool {
    if opening.index != index {
        return false;
    }
    let leaf = [opening.salt.as_slice(), &opening.value.to_bytes()].concat();
    verify_merkle_proof(&leaf, index, num_leaves, &opening.proof, root)
}

fn verify_opening_ext(
    opening: &MerkleOpeningExt,
    index: usize,
    num_leaves: usize,
    root: &[u8],
) -> bool {
    if opening.index != index {
        return false;
    }
    verify_merkle_proof(&opening.value.to_bytes(), index, num_leaves, &opening.proof, root)
}

/// Merkle root of unsalted extension-field leaves; matches `build_merkle_tree_ext`.
fn merkle_root_of_ext(values: &[Ext]) -> Vec<u8> {
    let leaves: Vec<Vec<u8>> = values.iter().map(|v| v.to_bytes().to_vec()).collect();
    MerkleTree::new(leaves).root().unwrap()
}

/// Derive the out-of-domain point in the extension field, matching the prover:
/// reject only the (negligible) base-field case. Bounded, so a malformed
/// transcript cannot spin forever.
fn derive_z_verifier(transcript: &mut FiatShamirTranscript) -> Option<Ext> {
    for _ in 0..64 {
        let z = transcript.squeeze_ext_challenge();
        if !z.is_base() {
            return Some(z);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::program::trace::ExecutionTrace;
    use crate::fibonacci::StarkProver;

    fn fibonacci_list(n: usize) -> Vec<u64> {
        let mut fibs = Vec::with_capacity(n);
        let mut a = 1u64;
        let mut b = 1u64;
        for _ in 0..n {
            fibs.push(a);
            let next = a.wrapping_add(b);
            a = b;
            b = next;
        }
        fibs
    }

    fn make_valid_proof() -> StarkProof {
        let mut trace = ExecutionTrace::new();
        let fib: Vec<BabyBear> = fibonacci_list(64).iter().map(|x| BabyBear::new(*x)).collect();
        trace.insert_column(fib);
        let prover = StarkProver::new(trace);
        prover.generate_proof(false)
    }

    #[test]
    fn test_verifier_accepts_valid_proof() {
        let proof = make_valid_proof();
        let verifier = StarkVerifier;
        assert!(verifier.verify(&proof), "Verifier should accept a valid proof");
    }

    #[test]
    fn test_masking_is_zero_knowledge() {
        // Two proofs of the same trace use fresh blinding, so their openings
        // (here the OOD evaluation t_z) differ, yet both verify.
        let p1 = make_valid_proof();
        let p2 = make_valid_proof();
        let verifier = StarkVerifier;
        assert!(verifier.verify(&p1) && verifier.verify(&p2));
        assert_ne!(p1.t_z, p2.t_z, "masking should randomize the openings");
    }

    #[test]
    fn test_verifier_rejects_bad_ood_value() {
        let mut proof = make_valid_proof();
        // Tamper with OOD trace evaluation → breaks constraint check C(z)=Q(z)*Z(z)
        proof.t_z = proof.t_z + Ext::one();
        let verifier = StarkVerifier;
        assert!(!verifier.verify(&proof), "Verifier should reject tampered OOD value");
    }

    #[test]
    fn test_verifier_rejects_bad_fri_final() {
        let mut proof = make_valid_proof();
        // Tamper with one final-layer value → breaks constancy + commitment + fold consistency
        proof.fri_final_layer[0] = proof.fri_final_layer[0] + Ext::one();
        let verifier = StarkVerifier;
        assert!(!verifier.verify(&proof), "Verifier should reject tampered FRI final layer");
    }

    #[test]
    fn test_verifier_rejects_bad_trace_commitment() {
        let mut proof = make_valid_proof();
        // Tamper with trace commitment → Merkle proofs fail
        proof.trace_commitment[0] ^= 0xff;
        let verifier = StarkVerifier;
        assert!(
            !verifier.verify(&proof),
            "Verifier should reject tampered trace commitment"
        );
    }

    #[test]
    fn test_verifier_rejects_bad_quotient_commitment() {
        let mut proof = make_valid_proof();
        // Tamper with quotient commitment → Merkle proofs fail AND
        // Fiat-Shamir transcript diverges
        proof.quotient_commitment[0] ^= 0xff;
        let verifier = StarkVerifier;
        assert!(
            !verifier.verify(&proof),
            "Verifier should reject tampered quotient commitment"
        );
    }

    #[test]
    fn test_verifier_rejects_bad_fri_commitment() {
        let mut proof = make_valid_proof();
        // Tamper with first FRI commitment → transcript diverges + Merkle fails
        proof.fri_commitments[0][0] ^= 0xff;
        let verifier = StarkVerifier;
        assert!(
            !verifier.verify(&proof),
            "Verifier should reject tampered FRI commitment"
        );
    }

    /// Malformed proofs must be rejected, not panic. Each of these reached an
    /// `assert!` or an out-of-bounds index before the parameter validation and
    /// the bounds checks were added.
    #[test]
    fn malformed_proofs_do_not_panic() {
        for mutate in [
            (|p: &mut StarkProof| p.trace_len = 0) as fn(&mut StarkProof),
            |p| p.trace_len = 3,
            |p| p.trace_len = 1 << 30,
            |p| p.lde_size = 12345,
            |p| p.fri_commitments.clear(),
            |p| p.fri_final_layer.clear(),
            |p| p.trace_commitment.clear(),
            |p| {
                for qp in p.query_proofs.iter_mut() {
                    qp.fri_openings.clear();
                }
            },
        ] {
            let mut proof = make_valid_proof();
            mutate(&mut proof);
            assert!(!StarkVerifier.verify(&proof));
        }
    }

    /// Openings are bound to the index the verifier derived, so a leaf cannot
    /// be moved to another query slot.
    #[test]
    fn index_substitution_is_rejected() {
        let mut proof = make_valid_proof();
        proof.query_proofs[0].trace_opening = proof.query_proofs[1].trace_opening.clone();
        assert!(!StarkVerifier.verify(&proof));
    }

    #[test]
    fn test_verifier_rejects_wrong_query_count() {
        let mut proof = make_valid_proof();
        // Remove a query proof → wrong number of queries
        proof.query_proofs.pop();
        let verifier = StarkVerifier;
        assert!(
            !verifier.verify(&proof),
            "Verifier should reject wrong number of query proofs"
        );
    }
}
