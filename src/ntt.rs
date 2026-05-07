// NTT for BabyBear field.
//
// This module exposes a CPU implementation that is always available and a CUDA
// implementation gated behind the `cuda` feature. The CUDA path is in a private
// inline submodule but its surface (`ntt_cuda`, `intt_cuda`, `cuda_available`,
// `CudaBuffer`) is re-exported flat at this module's level so callers can stay
// agnostic.

use crate::babybear::BabyBear;

// ── CPU NTT (always available) ──────────────────────────────────────────────

#[inline]
fn bit_reverse(mut x: usize, log_n: usize) -> usize {
    let mut result = 0;
    for _ in 0..log_n {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    result
}

/// In-place Cooley-Tukey NTT.
pub fn ntt(values: &mut [BabyBear], omega: BabyBear) {
    let n = values.len();
    assert!(n.is_power_of_two(), "NTT size must be power of 2");
    let log_n = n.trailing_zeros() as usize;

    for i in 0..n {
        let j = bit_reverse(i, log_n);
        if i < j {
            values.swap(i, j);
        }
    }

    let mut len = 2;
    while len <= n {
        let step = n / len;
        let w_len = omega.pow(step as u64);

        for i in (0..n).step_by(len) {
            let mut w = BabyBear::one();
            for j in 0..len / 2 {
                let u = values[i + j];
                let v = values[i + j + len / 2] * w;
                values[i + j] = u + v;
                values[i + j + len / 2] = u - v;
                w = w * w_len;
            }
        }
        len *= 2;
    }
}

/// In-place inverse NTT.
pub fn intt(values: &mut [BabyBear], omega: BabyBear) {
    let n = values.len();

    let inv_omega = omega.pow(n as u64 - 1);
    ntt(values, inv_omega);

    let inv_n = BabyBear::new(n as u64).inverse();
    for v in values.iter_mut() {
        *v = *v * inv_n;
    }
}

/// Generate roots of unity domain.
pub fn roots_of_unity_domain(n: usize) -> Vec<BabyBear> {
    assert!(n.is_power_of_two(), "Domain size must be power of 2");
    let log_n = n.trailing_zeros();
    let omega = BabyBear::get_root_of_unity(log_n);

    let mut domain = Vec::with_capacity(n);
    let mut cur = BabyBear::one();
    for _ in 0..n {
        domain.push(cur);
        cur = cur * omega;
    }
    domain
}

// ── CUDA NTT (feature = "cuda") ─────────────────────────────────────────────

#[cfg(feature = "cuda")]
mod cuda {
    use super::BabyBear;
    use std::ffi::CStr;
    use std::os::raw::c_char;

    type CudaError = i32;

    const CUDA_SUCCESS: CudaError = 0;

    #[link(name = "ntt_cuda", kind = "static")]
    unsafe extern "C" {
        fn cuda_copy_to_device(d_dest: *mut u64, h_src: *const u64, count: usize) -> CudaError;
        fn cuda_copy_from_device(h_dest: *mut u64, d_src: *const u64, count: usize) -> CudaError;
        fn cuda_malloc(d_ptr: *mut *mut u64, count: usize) -> CudaError;
        fn cuda_free(d_ptr: *mut u64) -> CudaError;
        fn cuda_get_error_string(error: CudaError) -> *const c_char;
        fn cudaGetDeviceCount(count: *mut i32) -> CudaError;

        // Persistent NTT context: caches twiddles + reusable device buffer per n.
        // Contexts are cached for the lifetime of the process (see `get_or_create_ctx`),
        // so `ntt_ctx_destroy` is only declared on the C side; Rust never calls it.
        fn ntt_ctx_create(n: u32) -> *mut std::ffi::c_void;
        fn ntt_run_inplace(ctx: *mut std::ffi::c_void, h_data: *mut u64);
        fn intt_run_inplace(ctx: *mut std::ffi::c_void, h_data: *mut u64);
    }

    // Compile-time guarantee: BabyBear has the same memory layout as u64 so we
    // can hand a `&mut [BabyBear]` straight to the C side as a `*mut u64` without
    // allocating an intermediate Vec.
    const _: () = assert!(std::mem::size_of::<BabyBear>() == std::mem::size_of::<u64>());
    const _: () = assert!(std::mem::align_of::<BabyBear>() == std::mem::align_of::<u64>());

    // Wrapper to give a raw NttCtx pointer Send/Sync. The pointer points at
    // device-managed state that the C side serialises internally; calls into it
    // from the prover are single-threaded today, but the cache itself is shared.
    struct CtxPtr(*mut std::ffi::c_void);
    unsafe impl Send for CtxPtr {}
    unsafe impl Sync for CtxPtr {}

    /// Look up (or lazily create) the cached NTT context for a given size.
    /// Subsequent calls for the same `n` skip device allocation, host twiddle
    /// computation, and per-stage H2D copies entirely.
    fn get_or_create_ctx(n: usize) -> *mut std::ffi::c_void {
        use std::collections::HashMap;
        use std::sync::{Mutex, OnceLock};

        static CACHE: OnceLock<Mutex<HashMap<usize, CtxPtr>>> = OnceLock::new();
        let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
        let mut guard = cache.lock().unwrap();
        if let Some(p) = guard.get(&n) {
            return p.0;
        }
        let ptr = unsafe { ntt_ctx_create(n as u32) };
        guard.insert(n, CtxPtr(ptr));
        ptr
    }

    /// Check if CUDA is available
    pub fn cuda_available() -> bool {
        unsafe {
            let mut count = 0;
            let err = cudaGetDeviceCount(&mut count);
            err == CUDA_SUCCESS && count > 0
        }
    }

    /// RAII wrapper for CUDA device memory
    pub struct CudaBuffer {
        ptr: *mut u64,
        size: usize,
    }

    impl CudaBuffer {
        pub fn new(size: usize) -> Result<Self, String> {
            let mut ptr: *mut u64 = std::ptr::null_mut();
            unsafe {
                let err = cuda_malloc(&mut ptr, size);
                if err != CUDA_SUCCESS {
                    let err_str = CStr::from_ptr(cuda_get_error_string(err));
                    return Err(format!("CUDA malloc failed: {}", err_str.to_string_lossy()));
                }
            }
            Ok(Self { ptr, size })
        }

        pub fn copy_from_host(&mut self, data: &[u64]) -> Result<(), String> {
            assert_eq!(data.len(), self.size, "Size mismatch");
            unsafe {
                let err = cuda_copy_to_device(self.ptr, data.as_ptr(), self.size);
                if err != CUDA_SUCCESS {
                    let err_str = CStr::from_ptr(cuda_get_error_string(err));
                    return Err(format!(
                        "CUDA copy to device failed: {}",
                        err_str.to_string_lossy()
                    ));
                }
            }
            Ok(())
        }

        pub fn copy_to_host(&self, data: &mut [u64]) -> Result<(), String> {
            assert_eq!(data.len(), self.size, "Size mismatch");
            unsafe {
                let err = cuda_copy_from_device(data.as_mut_ptr(), self.ptr, self.size);
                if err != CUDA_SUCCESS {
                    let err_str = CStr::from_ptr(cuda_get_error_string(err));
                    return Err(format!(
                        "CUDA copy from device failed: {}",
                        err_str.to_string_lossy()
                    ));
                }
            }
            Ok(())
        }

        pub fn as_ptr(&self) -> *mut u64 {
            self.ptr
        }
    }

    impl Drop for CudaBuffer {
        fn drop(&mut self) {
            unsafe {
                cuda_free(self.ptr);
            }
        }
    }

    unsafe impl Send for CudaBuffer {}
    unsafe impl Sync for CudaBuffer {}

    /// Perform NTT on GPU.
    ///
    /// Reinterprets `&mut [BabyBear]` as `*mut u64` (safe because BabyBear is
    /// `#[repr(C)] { value: u64 }`, statically asserted above) and hands it to a
    /// cached NTT context. That context owns: a reusable device buffer, forward
    /// twiddles, inverse twiddles, and n^-1 mod p — all reused across every call
    /// for a given size, so per-call overhead is one H2D + kernels + one D2H.
    pub fn ntt_cuda(values: &mut [BabyBear]) -> Result<(), String> {
        if !cuda_available() {
            return Err("CUDA not available".to_string());
        }
        let n = values.len();
        assert!(n.is_power_of_two(), "NTT size must be power of 2");
        assert!(n.trailing_zeros() <= 27, "BabyBear only supports NTT up to 2^27");

        let ctx = get_or_create_ctx(n);
        let raw = values.as_mut_ptr() as *mut u64;
        unsafe { ntt_run_inplace(ctx, raw); }
        Ok(())
    }

    /// Perform inverse NTT on GPU. See [`ntt_cuda`] for the caching model.
    pub fn intt_cuda(values: &mut [BabyBear]) -> Result<(), String> {
        if !cuda_available() {
            return Err("CUDA not available".to_string());
        }
        let n = values.len();
        assert!(n.is_power_of_two(), "NTT size must be power of 2");
        assert!(n.trailing_zeros() <= 27, "BabyBear only supports NTT up to 2^27");

        let ctx = get_or_create_ctx(n);
        let raw = values.as_mut_ptr() as *mut u64;
        unsafe { intt_run_inplace(ctx, raw); }
        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use super::super::{intt as cpu_intt, ntt as cpu_ntt};

        #[test]
        fn test_cuda_available() {
            println!("CUDA available: {}", cuda_available());
        }

        #[test]
        fn test_cuda_ntt_vs_cpu() {
            if !cuda_available() {
                println!("CUDA not available, skipping test");
                return;
            }

            let n = 256usize;
            let mut cpu_values: Vec<BabyBear> =
                (0..n).map(|i| BabyBear::new((i * 7 + 3) as u64)).collect();
            let mut gpu_values = cpu_values.clone();

            let omega = BabyBear::get_root_of_unity(n.trailing_zeros());

            cpu_ntt(&mut cpu_values, omega);
            ntt_cuda(&mut gpu_values).unwrap();

            for (i, (cpu_val, gpu_val)) in cpu_values.iter().zip(gpu_values.iter()).enumerate() {
                assert_eq!(
                    cpu_val.value, gpu_val.value,
                    "Mismatch at index {}: CPU={}, GPU={}",
                    i, cpu_val.value, gpu_val.value
                );
            }
        }

        #[test]
        fn test_cuda_intt_roundtrip() {
            if !cuda_available() {
                println!("CUDA not available, skipping test");
                return;
            }

            let n = 256usize;
            let original: Vec<BabyBear> = (0..n).map(|i| BabyBear::new((i * 7 + 3) as u64)).collect();
            let mut values = original.clone();

            ntt_cuda(&mut values).unwrap();
            intt_cuda(&mut values).unwrap();

            for (i, (orig, recovered)) in original.iter().zip(values.iter()).enumerate() {
                assert_eq!(
                    orig.value, recovered.value,
                    "Roundtrip failed at index {}",
                    i
                );
            }
        }
    }
}

#[cfg(feature = "cuda")]
pub use cuda::{cuda_available, intt_cuda, ntt_cuda, CudaBuffer};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_intt_roundtrip() {
        let n = 256usize;
        let omega = BabyBear::get_root_of_unity(n.trailing_zeros());

        let mut values: Vec<BabyBear> = (0..n).map(|i| BabyBear::new((i * 7 + 3) as u64)).collect();

        let original = values.clone();

        ntt(&mut values, omega);
        intt(&mut values, omega);

        for (a, b) in original.iter().zip(values.iter()) {
            assert_eq!(a.value, b.value, "NTT roundtrip failed");
        }
    }

    #[test]
    fn test_polynomial_evaluation() {
        let n = 8usize;
        let omega = BabyBear::get_root_of_unity(n.trailing_zeros());
        let domain = roots_of_unity_domain(n);

        let mut coeffs = vec![BabyBear::zero(); n];
        coeffs[0] = BabyBear::new(1);
        coeffs[1] = BabyBear::new(2);
        coeffs[2] = BabyBear::new(3);

        let mut evals = coeffs.clone();
        ntt(&mut evals, omega);

        assert_eq!(evals[0].value, 6);

        let x = domain[1];
        let expected = coeffs[0] + coeffs[1] * x + coeffs[2] * x * x;
        assert_eq!(evals[1].value, expected.value);
    }

    #[test]
    fn test_roots_of_unity() {
        let n = 16usize;
        let domain = roots_of_unity_domain(n);

        assert_eq!(domain[0].value, 1);

        let omega = domain[1];
        let result = omega.pow(n as u64);
        assert_eq!(result.value, 1);

        for i in 0..n {
            for j in i + 1..n {
                assert_ne!(
                    domain[i].value, domain[j].value,
                    "domain[{}] = {} equals domain[{}] = {}",
                    i, domain[i].value, j, domain[j].value
                );
            }
        }
    }
}
