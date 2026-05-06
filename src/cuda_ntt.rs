// CUDA NTT FFI bindings and safe Rust wrapper
// Ported from joda

use crate::babybear::BabyBear;
use std::ffi::CStr;
use std::os::raw::c_char;

type CudaError = i32;

const CUDA_SUCCESS: CudaError = 0;

#[link(name = "ntt_cuda", kind = "static")]
unsafe extern "C" {
    fn cuda_ntt(d_values: *mut u64, n: u32, omega: u64);
    fn cuda_intt(d_values: *mut u64, n: u32, omega: u64);
    fn cuda_copy_to_device(d_dest: *mut u64, h_src: *const u64, count: usize) -> CudaError;
    fn cuda_copy_from_device(h_dest: *mut u64, d_src: *const u64, count: usize) -> CudaError;
    fn cuda_malloc(d_ptr: *mut *mut u64, count: usize) -> CudaError;
    fn cuda_free(d_ptr: *mut u64) -> CudaError;
    fn cuda_get_error_string(error: CudaError) -> *const c_char;
    fn cudaGetDeviceCount(count: *mut i32) -> CudaError;

    pub fn cuda_ntt_batched(
        d_values: *mut u64,
        num_ntts: u32,
        ntt_size: u32,
        stride: u32,
        omega: u64,
    );
    pub fn cuda_intt_batched(
        d_values: *mut u64,
        num_ntts: u32,
        ntt_size: u32,
        stride: u32,
        omega: u64,
    );

    pub fn cuda_rs_encode_vertical(
        d_input: *const u64,
        d_output: *mut u64,
        d_intt_work: *mut u64,
        num_positions: u32,
        ntt_size_k: u32,
        ntt_size_kn: u32,
        omega_k: u64,
        omega_kn: u64,
    );

    // Persistent NTT context: caches twiddles + reusable device buffer per n.
    fn ntt_ctx_create(n: u32) -> *mut std::ffi::c_void;
    fn ntt_ctx_destroy(ctx: *mut std::ffi::c_void);
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
    use crate::ntt_babybear::{intt as cpu_intt, ntt as cpu_ntt};

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
