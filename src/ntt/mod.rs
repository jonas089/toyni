//! Number-theoretic transform.
//!
//! `cpu` is the reference implementation and the source of truth for
//! correctness. `cuda` wraps the GPU kernel in `cuda/ntt_kernel.cu` and is
//! only compiled with the `cuda` feature.

pub mod cpu;

#[cfg(feature = "cuda")]
pub mod cuda;

pub use cpu::{intt, intt_babybear, ntt, ntt_babybear, roots_of_unity};

#[cfg(feature = "cuda")]
pub use cuda::{cuda_available, intt_cuda, ntt_cuda, PinnedSquare};
