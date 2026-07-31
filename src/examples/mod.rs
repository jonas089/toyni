//! Example AIRs, generic over the base field so each runs on both engines.
//!
//! - [`fibonacci`]: single-column second-order recurrence (the baseline).
//! - [`range_check`]: multi-column byte decomposition + running sum.
//! - [`hash_chain`]: Poseidon-style x^5 permutation chain (heavy workload).

pub mod fibonacci;
pub mod hash_chain;
pub mod range_check;

pub use fibonacci::FibonacciAir;
pub use hash_chain::HashChainAir;
pub use range_check::RangeCheckAir;
