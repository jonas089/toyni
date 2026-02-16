use crate::babybear::BabyBear;
use crate::math::domain::BabyBearDomain;
use crate::math::polynomial::Polynomial;
use crate::merkle::verify_merkle_proof;
use crate::prover::{
    eval_boundary_1, eval_boundary_2, eval_fibonacci_constraint, MerkleOpening, StarkProof,
    BLOWUP, COSET_SHIFT, NUM_QUERIES,
};
use crate::transcript::FiatShamirTranscript;

pub struct StarkVerifier;

impl StarkVerifier {
    pub fn verify(&self, proof: &StarkProof) -> bool {
        let trace_len = proof.trace_len;
        let lde_size = proof.lde_size;

        // Sanity: lde_size must equal trace_len * BLOWUP
        if lde_size != trace_len * BLOWUP {
            return false;
        }

        let domain = BabyBearDomain::new(trace_len);
        let extended_domain = BabyBearDomain::new(lde_size);
        let shift = BabyBear::new(COSET_SHIFT);
        let shifted_domain = extended_domain.get_coset(shift);
        let g = domain.group_gen();

        let z_poly = Polynomial::new(domain.vanishing_poly_coeffs());

        // ── 1. Replay Fiat-Shamir transcript ───────────────────────────
        let mut transcript = FiatShamirTranscript::new();
        transcript.absorb_commitment(&proof.trace_commitment);
        transcript.absorb_commitment(&proof.quotient_commitment);

        let z = derive_z_verifier(&mut transcript, &extended_domain, &shifted_domain);

        transcript.absorb_field(proof.t_z);
        transcript.absorb_field(proof.t_gz);
        transcript.absorb_field(proof.t_ggz);
        transcript.absorb_field(proof.q_z);

        // ── 2. OOD constraint check: C(z) = Q(z) · Z(z) ──────────────
        let c_z = eval_fibonacci_constraint(proof.t_ggz, proof.t_gz, proof.t_z)
            * eval_boundary_1(z, g, trace_len)
            * eval_boundary_2(z, g, trace_len);
        if c_z != proof.q_z * z_poly.evaluate(z) {
            return false;
        }

        // ── 3. Replay FRI commitments & derive betas ───────────────────
        if proof.fri_commitments.is_empty() {
            return false;
        }
        transcript.absorb_commitment(&proof.fri_commitments[0]);

        let num_fri_folds = proof.fri_commitments.len() - 1;
        let mut fri_betas = Vec::with_capacity(num_fri_folds);

        for i in 1..proof.fri_commitments.len() {
            let beta = transcript.squeeze_challenge();
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

            // 6a. Verify Merkle proofs for trace openings (3 positions)
            if !verify_opening(&qp.trace_opening, &proof.trace_commitment) {
                return false;
            }
            if !verify_opening(&qp.trace_opening_g, &proof.trace_commitment) {
                return false;
            }
            if !verify_opening(&qp.trace_opening_gg, &proof.trace_commitment) {
                return false;
            }

            // Verify indices are correct
            let expected_idx_g = (qi + BLOWUP) % lde_size;
            let expected_idx_gg = (qi + 2 * BLOWUP) % lde_size;
            if qp.trace_opening.index != qi
                || qp.trace_opening_g.index != expected_idx_g
                || qp.trace_opening_gg.index != expected_idx_gg
            {
                return false;
            }

            // 6b. Verify Merkle proof for quotient opening
            if !verify_opening(&qp.quotient_opening, &proof.quotient_commitment) {
                return false;
            }

            // 6c. Verify Merkle proofs for DEEP layer (FRI layer 0)
            if !verify_opening(&qp.deep_opening, &proof.fri_commitments[0]) {
                return false;
            }
            if !verify_opening(&qp.deep_opening_pair, &proof.fri_commitments[0]) {
                return false;
            }

            // 6d. DEEP polynomial consistency check
            //     D(x) = (Q(x)-Q(z))/(x-z) + (T(x)-T(z))/(x-z)
            //           + (T(gx)-T(gz))/(x-z) + (T(g²x)-T(g²z))/(x-z)
            //
            //     Verify the opened DEEP value matches the reconstruction
            //     from the opened trace and quotient values + OOD values.
            let x_i = shifted_elements[qi];
            let t_x = qp.trace_opening.value;
            let t_gx = qp.trace_opening_g.value;
            let t_ggx = qp.trace_opening_gg.value;
            let q_x = qp.quotient_opening.value;

            let inv_x_minus_z = (x_i - z).inverse();
            let expected_deep = (q_x - proof.q_z) * inv_x_minus_z
                + (t_ggx - proof.t_ggz) * inv_x_minus_z
                + (t_gx - proof.t_gz) * inv_x_minus_z
                + (t_x - proof.t_z) * inv_x_minus_z;

            if qp.deep_opening.value != expected_deep {
                return false;
            }

            // 6e. First FRI fold: layer 0 → layer 1
            let a0 = qp.deep_opening.value;
            let b0 = qp.deep_opening_pair.value;
            let x0 = shifted_elements[qi];

            let mut prev_folded = {
                let avg = (a0 + b0) * half_inv;
                let diff = (a0 - b0) * half_inv;
                avg + diff * fri_betas[0] * x0.inverse()
            };

            // 6f. Intermediate FRI layers
            let mut pos = qi;

            for layer in 0..qp.fri_openings.len() {
                let fold_k = layer + 1;
                let layer_size = lde_size >> fold_k;
                let half = layer_size / 2;

                let lo = pos % half;
                let in_first_half = pos == lo;

                let (ref op, ref op_pair) = qp.fri_openings[layer];

                // Merkle proofs
                if !verify_opening(op, &proof.fri_commitments[fold_k]) {
                    return false;
                }
                if !verify_opening(op_pair, &proof.fri_commitments[fold_k]) {
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
                let x = shifted_elements[lo].pow(1u64 << fold_k);

                let a_l = op.value;
                let b_l = op_pair.value;
                let avg = (a_l + b_l) * half_inv;
                let diff = (a_l - b_l) * half_inv;
                prev_folded = avg + diff * fri_betas[fold_k] * x.inverse();

                pos = lo;
            }

            // 6g. Final FRI value check
            if prev_folded != proof.fri_final_value {
                return false;
            }
        }

        true
    }
}

fn verify_opening(opening: &MerkleOpening, root: &[u8]) -> bool {
    let leaf = opening.value.to_bytes().to_vec();
    verify_merkle_proof(leaf, &opening.proof, &root.to_vec())
}

fn derive_z_verifier(
    transcript: &mut FiatShamirTranscript,
    extended_domain: &BabyBearDomain,
    shifted_domain: &BabyBearDomain,
) -> BabyBear {
    let ext_set: std::collections::HashSet<BabyBear> =
        extended_domain.elements().into_iter().collect();
    let shift_set: std::collections::HashSet<BabyBear> =
        shifted_domain.elements().into_iter().collect();
    let g = extended_domain.group_gen();

    loop {
        let z = transcript.squeeze_challenge();
        if !ext_set.contains(&z)
            && !shift_set.contains(&z)
            && !shift_set.contains(&(g * z))
            && !shift_set.contains(&(g * g * z))
        {
            return z;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::program::trace::ExecutionTrace;
    use crate::prover::StarkProver;

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
    fn test_verifier_rejects_bad_ood_value() {
        let mut proof = make_valid_proof();
        // Tamper with OOD trace evaluation → breaks constraint check C(z)=Q(z)*Z(z)
        proof.t_z = proof.t_z + BabyBear::one();
        let verifier = StarkVerifier;
        assert!(!verifier.verify(&proof), "Verifier should reject tampered OOD value");
    }

    #[test]
    fn test_verifier_rejects_bad_fri_final() {
        let mut proof = make_valid_proof();
        // Tamper with final FRI constant → breaks fold consistency
        proof.fri_final_value = proof.fri_final_value + BabyBear::one();
        let verifier = StarkVerifier;
        assert!(!verifier.verify(&proof), "Verifier should reject tampered FRI final value");
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
