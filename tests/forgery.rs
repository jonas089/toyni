//! END-TO-END FORGERY PoC.
//!
//! Claim under test: the fold-to-1 FRI in `fibonacci.rs` / `verifier.rs` does
//! not enforce ANY low-degree bound, and this is exploitable end-to-end — a
//! malicious prover can make `StarkVerifier::verify` accept a proof of a trace
//! that VIOLATES the Fibonacci recurrence.
//!
//! The exploit:
//!   1. Take a trace that breaks the recurrence. Then C(x) = constraint·boundary
//!      does NOT vanish on H, so the quotient q(x) = C(x)/Z(x) is a rational
//!      function — as a vector on the LDE coset it is FAR from any low-degree
//!      polynomial. The honest prover panics on its own OOD self-check; a
//!      malicious prover simply omits that check.
//!   2. Forge q_z := C(z)/Z(z) so the verifier's single-point OOD check
//!      C(z) == q_z·Z(z) holds by construction.
//!   3. Build the DEEP layer exactly as the verifier reconstructs it, so the
//!      per-query DEEP consistency check passes.
//!   4. Fold the (high-degree) DEEP codeword honestly all the way to length 1.
//!      Every fold-consistency check + the final-value check pass, because the
//!      prover folded honestly — fold-consistency never tests degree.
//!
//! BEFORE the fix (fold-to-1 FRI, single scalar final check) this forgery was
//! accepted by `StarkVerifier::verify` — a confirmed end-to-end soundness break.
//! AFTER the fix (fold to the size-BLOWUP degree-bound layer + final-layer
//! constancy check) the high-degree DEEP codeword folds to a NON-constant final
//! layer and is rejected. This test now serves as a regression guard: it asserts
//! the forgery is REJECTED.

use toyni::babybear::BabyBear;
use toyni::fibonacci::{
    eval_boundary_1, eval_boundary_2, eval_fibonacci_constraint, MerkleOpening, QueryProof,
    StarkProof, StarkProver, BLOWUP, COSET_SHIFT, MASK_DEGREE, NUM_QUERIES,
};
use toyni::math::domain::BabyBearDomain;
use toyni::math::fri::fri_fold;
use toyni::math::polynomial::Polynomial;
use toyni::merkle::MerkleTree;
use toyni::program::trace::ExecutionTrace;
use toyni::transcript::FiatShamirTranscript;
use toyni::verifier::StarkVerifier;

fn build_merkle_tree(evals: &[BabyBear]) -> MerkleTree {
    let leaves: Vec<Vec<u8>> = evals.iter().map(|v| v.to_bytes().to_vec()).collect();
    MerkleTree::new(leaves)
}

fn open_merkle(tree: &MerkleTree, evals: &[BabyBear], index: usize) -> MerkleOpening {
    let proof = tree.get_proof(index).expect("Index out of bounds");
    MerkleOpening {
        index,
        value: evals[index],
        proof,
    }
}

// Replicates `derive_z_from_transcript` / `derive_z_verifier` verbatim so the
// forged transcript stays in lock-step with the verifier's replay.
fn derive_z(
    transcript: &mut FiatShamirTranscript,
    extended_domain: &BabyBearDomain,
    shifted_domain: &BabyBearDomain,
) -> BabyBear {
    use std::collections::HashSet;
    let ext_set: HashSet<BabyBear> = extended_domain.elements().into_iter().collect();
    let shift_set: HashSet<BabyBear> = shifted_domain.elements().into_iter().collect();
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

/// Malicious prover: same pipeline as `StarkProver::generate_proof` EXCEPT
///   - no OOD self-check assertion, and
///   - q_z is forged as C(z)/Z(z) instead of q_poly(z).
fn forge_proof(trace_col: Vec<BabyBear>) -> StarkProof {
    let trace_len = trace_col.len();
    let domain = BabyBearDomain::new(trace_len);
    let lde_size = trace_len * BLOWUP;
    let extended_domain = BabyBearDomain::new(lde_size);
    let shift = BabyBear::new(COSET_SHIFT);
    let shifted_domain = extended_domain.get_coset(shift);

    let z_poly = Polynomial::new(domain.vanishing_poly_coeffs());
    let g = domain.group_gen();

    // 1. trace polynomial (honest interpolation of the *invalid* trace)
    let domain_elements = domain.elements();
    let mut et = ExecutionTrace::new();
    et.insert_column(trace_col);
    let trace_poly = et.interpolate_column(&domain_elements, 0);

    let shifted_elements = shifted_domain.elements();
    let trace_lde: Vec<BabyBear> = shifted_elements
        .iter()
        .map(|&x| trace_poly.evaluate(x))
        .collect();
    let trace_tree = build_merkle_tree(&trace_lde);
    let trace_commitment = trace_tree.root().unwrap();

    // 2. constraint & quotient. For an invalid trace C does NOT vanish on H,
    //    so q_evals is a rational function's values: high-degree as a vector.
    let c_evals: Vec<BabyBear> = shifted_elements
        .iter()
        .map(|&x| {
            eval_fibonacci_constraint(
                trace_poly.evaluate(g * g * x),
                trace_poly.evaluate(g * x),
                trace_poly.evaluate(x),
            ) * eval_boundary_1(x, g, trace_len)
                * eval_boundary_2(x, g, trace_len)
        })
        .collect();
    let q_evals: Vec<BabyBear> = shifted_elements
        .iter()
        .zip(c_evals.iter())
        .map(|(&x, &c)| c / z_poly.evaluate(x))
        .collect();
    let quotient_tree = build_merkle_tree(&q_evals);
    let quotient_commitment = quotient_tree.root().unwrap();

    // 3. Fiat-Shamir → z
    let mut transcript = FiatShamirTranscript::new();
    transcript.absorb_commitment(&trace_commitment);
    transcript.absorb_commitment(&quotient_commitment);
    let z = derive_z(&mut transcript, &extended_domain, &shifted_domain);

    // 4. OOD evaluations. t_* are honest; q_z is FORGED so C(z) == q_z·Z(z).
    let t_z = trace_poly.evaluate(z);
    let t_gz = trace_poly.evaluate(g * z);
    let t_ggz = trace_poly.evaluate(g * g * z);
    let c_z = eval_fibonacci_constraint(t_ggz, t_gz, t_z)
        * eval_boundary_1(z, g, trace_len)
        * eval_boundary_2(z, g, trace_len);
    let q_z = c_z / z_poly.evaluate(z); // <-- the forgery

    transcript.absorb_field(t_z);
    transcript.absorb_field(t_gz);
    transcript.absorb_field(t_ggz);
    transcript.absorb_field(q_z);

    // 5. DEEP layer, reconstructed exactly as the verifier does (index-based),
    //    using the forged q_z so the DEEP consistency check matches.
    let d_evals: Vec<BabyBear> = (0..lde_size)
        .map(|i| {
            let x = shifted_elements[i];
            let t_x = trace_lde[i];
            let t_gx = trace_lde[(i + BLOWUP) % lde_size];
            let t_ggx = trace_lde[(i + 2 * BLOWUP) % lde_size];
            let q_x = q_evals[i];
            let ixz = (x - z).inverse();
            (q_x - q_z) * ixz + (t_ggx - t_ggz) * ixz + (t_gx - t_gz) * ixz + (t_x - t_z) * ixz
        })
        .collect();

    // 6. FRI fold to length 1 (verbatim from the prover)
    let mut fri_layers: Vec<Vec<BabyBear>> = Vec::new();
    let mut fri_trees: Vec<MerkleTree> = Vec::new();
    let mut fri_commitments: Vec<Vec<u8>> = Vec::new();
    fri_layers.push(d_evals.clone());
    let tree0 = build_merkle_tree(&d_evals);
    let root0 = tree0.root().unwrap();
    transcript.absorb_commitment(&root0);
    fri_commitments.push(root0);
    fri_trees.push(tree0);

    // The malicious prover follows the FIXED protocol as closely as it can:
    // fold honestly down to the degree-bound layer. Because the DEEP codeword
    // is high-degree (invalid trace ⇒ rational quotient), this final layer is
    // NOT constant, so the verifier's final-layer low-degree check rejects it.
    let final_layer_size = lde_size / (trace_len + MASK_DEGREE).next_power_of_two();
    let mut current = d_evals;
    let mut xs: Vec<BabyBear> = shifted_elements.clone();
    while current.len() > final_layer_size {
        let beta = transcript.squeeze_challenge();
        let folded = fri_fold(&current, &xs, beta);
        xs.truncate(folded.len());
        for x in &mut xs {
            *x = *x * *x;
        }
        fri_layers.push(folded.clone());
        let tree = build_merkle_tree(&folded);
        let root = tree.root().unwrap();
        transcript.absorb_commitment(&root);
        fri_commitments.push(root);
        fri_trees.push(tree);
        current = folded;
    }
    let fri_final_layer = current;

    // 7. query phase (verbatim from the prover)
    let first_layer_half = fri_layers[0].len() / 2;
    let query_indices = transcript.squeeze_indices(NUM_QUERIES, first_layer_half);
    let mut query_proofs = Vec::with_capacity(NUM_QUERIES);
    for &qi in &query_indices {
        let idx_g = (qi + BLOWUP) % lde_size;
        let idx_gg = (qi + 2 * BLOWUP) % lde_size;
        let trace_opening = open_merkle(&trace_tree, &trace_lde, qi);
        let trace_opening_g = open_merkle(&trace_tree, &trace_lde, idx_g);
        let trace_opening_gg = open_merkle(&trace_tree, &trace_lde, idx_gg);
        let quotient_opening = open_merkle(&quotient_tree, &q_evals, qi);

        let half0 = fri_layers[0].len() / 2;
        let deep_opening = open_merkle(&fri_trees[0], &fri_layers[0], qi);
        let deep_opening_pair = open_merkle(&fri_trees[0], &fri_layers[0], qi + half0);

        let mut fri_openings = Vec::new();
        let mut idx = qi;
        for layer_idx in 1..fri_layers.len() - 1 {
            let half = fri_layers[layer_idx].len() / 2;
            idx %= half;
            let op = open_merkle(&fri_trees[layer_idx], &fri_layers[layer_idx], idx);
            let op_pair = open_merkle(&fri_trees[layer_idx], &fri_layers[layer_idx], idx + half);
            fri_openings.push((op, op_pair));
        }

        query_proofs.push(QueryProof {
            index: qi,
            deep_opening,
            deep_opening_pair,
            trace_opening,
            trace_opening_g,
            trace_opening_gg,
            quotient_opening,
            fri_openings,
        });
    }

    StarkProof {
        trace_len,
        lde_size,
        trace_commitment,
        quotient_commitment,
        t_z,
        t_gz,
        t_ggz,
        q_z,
        fri_commitments,
        fri_final_layer,
        query_proofs,
    }
}

fn fib_list(n: usize) -> Vec<u64> {
    let mut fibs = Vec::with_capacity(n);
    let (mut a, mut b) = (1u64, 1u64);
    for _ in 0..n {
        fibs.push(a);
        let next = a.wrapping_add(b);
        a = b;
        b = next;
    }
    fibs
}

#[test]
fn forged_invalid_trace_is_rejected_by_fixed_fri() {
    // Start from a valid Fibonacci trace, then corrupt one interior value so the
    // recurrence is violated at rows ~28-30 (well inside the enforced range).
    let mut fib = fib_list(64);
    fib[30] = fib[30].wrapping_add(12345);
    let trace: Vec<BabyBear> = fib.iter().map(|&v| BabyBear::new(v)).collect();

    // (a) The HONEST prover refuses this trace: its OOD self-check panics.
    let prev_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let honest = std::panic::catch_unwind(|| {
        let mut et = ExecutionTrace::new();
        et.insert_column(trace.clone());
        StarkProver::new(et).generate_proof(false)
    });
    std::panic::set_hook(prev_hook);
    assert!(
        honest.is_err(),
        "expected the honest prover to reject the invalid trace (it did not — the trace may be valid)"
    );
    println!("honest prover rejected the invalid trace: {}", honest.is_err());

    // (b) The MALICIOUS prover forges a proof of the SAME invalid trace.
    let forged = forge_proof(trace);
    let accepted = StarkVerifier.verify(&forged);
    println!("forged proof of an INVALID trace accepted by FIXED verifier: {accepted}");

    assert!(
        !accepted,
        "REGRESSION: the final-layer low-degree FRI check must REJECT this forgery \
         (a high-degree DEEP codeword folds to a NON-constant size-BLOWUP final layer)"
    );
}
