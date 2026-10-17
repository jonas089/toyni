//! TOYNI-04 regression: trace_len = 4 is rejected by the TOYNI-01 rate check before any division.

use toyni::ext::Ext;
use toyni::fibonacci::StarkProof;
use toyni::verifier::StarkVerifier;

#[test]
fn verifier_rejects_trace_len_4_without_panicking() {
    // All-zero OOD claims satisfy C(z) = q(z) Z_H(z) trivially (0 = 0).
    let proof = StarkProof {
        trace_len: 4,
        lde_size: 128,
        trace_commitment: vec![0; 32],
        quotient_commitment: vec![0; 32],
        t_z: Ext::zero(),
        t_gz: Ext::zero(),
        t_ggz: Ext::zero(),
        q_z: Ext::zero(),
        fri_commitments: vec![vec![0; 32]],
        fri_final_layer: vec![],
        query_proofs: vec![],
    };
    let result = std::panic::catch_unwind(|| StarkVerifier.verify(&proof));
    assert_eq!(result.ok(), Some(false));
}
