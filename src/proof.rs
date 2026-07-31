//! Proof data structures and protocol parameters, generic over the engine.

use crate::engine::Engine;
use crate::merkle::Hash;

/// Protocol parameters. Defaults give a tested Reed-Solomon rate of `2^-3`
/// with 44 queries plus ~124-bit challenge-field soundness.
#[derive(Debug, Clone, Copy)]
pub struct ProofOptions {
    /// log2 of the inverse rate above the FRI degree-bound headroom.
    pub log_blowup: u32,
    /// Number of FRI spot-check queries.
    pub num_queries: usize,
    /// Zero-knowledge: blind every committed column and salt all leaves.
    pub zk: bool,
}

impl Default for ProofOptions {
    fn default() -> Self {
        Self {
            log_blowup: 3,
            num_queries: 44,
            zk: true,
        }
    }
}

impl ProofOptions {
    /// log2 of the number of masking coefficients per column: enough to cover
    /// every revealed evaluation (2 per query + the OOD mask points), capped so
    /// the masked degree stays within the FRI bound.
    pub fn mask_log_size(&self, log_trace: u32) -> u32 {
        let want = (2 * self.num_queries + 4).next_power_of_two().trailing_zeros();
        want.min(log_trace.saturating_sub(2))
    }
}

/// Opening of one committed trace position: the column values, salt, path.
#[derive(Debug, Clone)]
pub struct TraceOpening<E: Engine> {
    pub values: Vec<E::Base>,
    pub salt: Vec<u8>,
    pub path: Vec<Hash>,
}

/// Opening of one committed composition position.
#[derive(Debug, Clone)]
pub struct CompositionOpening<E: Engine> {
    pub value: E::Ext,
    pub salt: Vec<u8>,
    pub path: Vec<Hash>,
}

/// All openings for one query position `q` and its sibling `q + |D|/2`.
#[derive(Debug, Clone)]
pub struct QueryOpenings<E: Engine> {
    pub trace: TraceOpening<E>,
    pub trace_sibling: TraceOpening<E>,
    pub composition: CompositionOpening<E>,
    pub composition_sibling: CompositionOpening<E>,
}

/// A complete STARK proof.
#[derive(Debug, Clone)]
pub struct StarkProof<E: Engine> {
    pub log_trace_len: u32,
    pub trace_root: Hash,
    pub composition_root: Hash,
    /// `ood_trace[k][c]` = column `c` at the mask point `mask_point(gamma, k)`.
    pub ood_trace: Vec<Vec<E::Ext>>,
    /// The out-of-domain composition value `F(gamma)`.
    pub ood_composition: E::Ext,
    pub fri: E::FriProof,
    pub queries: Vec<QueryOpenings<E>>,
}
