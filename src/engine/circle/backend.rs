//! Compute-backend dispatch: every heavy prover stage funnels through these
//! wrappers, which route to the Metal GPU when the `metal` feature is enabled,
//! the device is available, and the problem is large enough — and to the
//! (rayon-parallel) CPU implementation otherwise. Results are bit-identical
//! either way.

use crate::engine::circle::cfft::{self, Twiddles};
use crate::field::{M31, QM31};
use crate::merkle::MerkleTree;

/// Forward CFFT (coefficients → values), in place.
pub fn evaluate_m31(values: &mut [M31], twiddles: &Twiddles) {
    #[cfg(feature = "metal")]
    if values.len() >= crate::engine::circle::metal::GPU_THRESHOLD
        && crate::engine::circle::metal::with_context(|ctx| ctx.evaluate_m31(values, twiddles)).is_some()
    {
        return;
    }
    cfft::evaluate(values, twiddles);
}

/// Inverse CFFT (values → coefficients), in place.
pub fn interpolate_m31(values: &mut [M31], twiddles: &Twiddles) {
    #[cfg(feature = "metal")]
    if values.len() >= crate::engine::circle::metal::GPU_THRESHOLD
        && crate::engine::circle::metal::with_context(|ctx| ctx.interpolate_m31(values, twiddles)).is_some()
    {
        return;
    }
    cfft::interpolate(values, twiddles);
}

/// Inverse CFFT over QM31 values, in place.
pub fn interpolate_qm31(values: &mut [QM31], twiddles: &Twiddles) {
    #[cfg(feature = "metal")]
    if values.len() >= crate::engine::circle::metal::GPU_THRESHOLD
        && crate::engine::circle::metal::with_context(|ctx| ctx.interpolate_qm31(values, twiddles)).is_some()
    {
        return;
    }
    cfft::interpolate(values, twiddles);
}

/// Zero-knowledge mask application: `lde[i] += vn[i] · r[i]`.
pub fn mask_add(lde: &mut [M31], vn: &[M31], r: &[M31]) {
    #[cfg(feature = "metal")]
    if lde.len() >= crate::engine::circle::metal::GPU_THRESHOLD
        && crate::engine::circle::metal::with_context(|ctx| ctx.mask_add(lde, vn, r)).is_some()
    {
        return;
    }
    for i in 0..lde.len() {
        lde[i] += vn[i] * r[i];
    }
}

/// Merkle tree over packed fixed-size leaves.
pub fn merkle_packed(packed: &[u8], leaf_len: usize) -> MerkleTree {
    let count = packed.len() / leaf_len;
    #[cfg(feature = "metal")]
    if count >= crate::engine::circle::metal::GPU_THRESHOLD && leaf_len < 159 {
        if let Some(tree) =
            crate::engine::circle::metal::with_context(|ctx| ctx.merkle(packed, leaf_len, count))
        {
            return tree;
        }
    }
    let _ = count;
    MerkleTree::from_packed(packed, leaf_len)
}

/// Merkle tree over QM31 evaluations (16-byte leaves, unsalted).
pub fn merkle_qm31(values: &[QM31]) -> MerkleTree {
    let mut packed = Vec::with_capacity(values.len() * 16);
    for v in values {
        packed.extend_from_slice(&v.to_bytes());
    }
    merkle_packed(&packed, 16)
}

/// One FRI fold layer.
pub fn fri_fold(values: &[QM31], inv_twiddles: &[M31], lambda: QM31) -> Vec<QM31> {
    #[cfg(feature = "metal")]
    if values.len() >= crate::engine::circle::metal::GPU_THRESHOLD {
        if let Some(out) =
            crate::engine::circle::metal::with_context(|ctx| ctx.fri_fold(values, inv_twiddles, lambda))
        {
            return out;
        }
    }
    crate::engine::circle::fri::fold_layer_cpu(values, inv_twiddles, lambda)
}
