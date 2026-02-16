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

/// Perform NTT on GPU
pub fn ntt_cuda(values: &mut [BabyBear]) -> Result<(), String> {
    if !cuda_available() {
        return Err("CUDA not available".to_string());
    }

    let n = values.len();
    assert!(n.is_power_of_two(), "NTT size must be power of 2");

    let omega = BabyBear::get_root_of_unity(n.trailing_zeros());

    let mut raw_values: Vec<u64> = values.iter().map(|v| v.value).collect();

    let mut d_buffer = CudaBuffer::new(n)?;
    d_buffer.copy_from_host(&raw_values)?;

    unsafe {
        cuda_ntt(d_buffer.as_ptr(), n as u32, omega.value);
    }

    d_buffer.copy_to_host(&mut raw_values)?;

    for (i, val) in raw_values.iter().enumerate() {
        values[i] = BabyBear::new(*val);
    }

    Ok(())
}

/// Perform inverse NTT on GPU
pub fn intt_cuda(values: &mut [BabyBear]) -> Result<(), String> {
    if !cuda_available() {
        return Err("CUDA not available".to_string());
    }

    let n = values.len();
    assert!(n.is_power_of_two(), "NTT size must be power of 2");

    let omega = BabyBear::get_root_of_unity(n.trailing_zeros());

    let mut raw_values: Vec<u64> = values.iter().map(|v| v.value).collect();

    let mut d_buffer = CudaBuffer::new(n)?;
    d_buffer.copy_from_host(&raw_values)?;

    unsafe {
        cuda_intt(d_buffer.as_ptr(), n as u32, omega.value);
    }

    d_buffer.copy_to_host(&mut raw_values)?;

    for (i, val) in raw_values.iter().enumerate() {
        values[i] = BabyBear::new(*val);
    }

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
