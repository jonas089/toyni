// CUDA NTT: FFI bindings and safe Rust wrapper.
//
// The C side (cuda/ntt_kernel.cu) owns all device memory. Values are canonical
// u32 BabyBear, one polynomial after another: polynomial c is the `n` values
// starting at `values[c * stride]`. Use `PinnedSquare` for the buffer so the
// copies overlap with the kernels and nothing is converted on the way in or out.

use crate::babybear::{BabyBear, BABYBEAR_PRIME};
use crate::ntt::roots_of_unity;
use std::ops::{Deref, DerefMut};

type CudaError = i32;
const CUDA_SUCCESS: CudaError = 0;

#[link(name = "ntt_cuda", kind = "static")]
unsafe extern "C" {
    fn cuda_ntt(values: *mut u32, roots: *const u32, n: u32, polys: u32, stride: u32) -> CudaError;
    fn cuda_intt(values: *mut u32, roots: *const u32, n: u32, polys: u32, stride: u32) -> CudaError;
    fn cuda_host_alloc(count: usize) -> *mut u32;
    fn cuda_host_free(p: *mut u32);
    fn cudaGetDeviceCount(count: *mut i32) -> CudaError;
}

/// Check if a CUDA device is present.
pub fn cuda_available() -> bool {
    let mut count = 0;
    let err = unsafe { cudaGetDeviceCount(&mut count) };
    err == CUDA_SUCCESS && count > 0
}

/// `cols` polynomials of `rows` canonical u32 field elements each, laid end to end
/// in pinned host memory. Column `c` (one ZODA column, one polynomial) is the slice
/// `self[c * rows .. (c + 1) * rows]`, so element `(row, col)` is at `col * rows + row`.
pub struct PinnedSquare {
    ptr: *mut u32,
    rows: usize,
    cols: usize,
}

impl PinnedSquare {
    /// Zero-filled square in pinned memory.
    pub fn new(rows: usize, cols: usize) -> Result<Self, String> {
        let count = (rows * cols).max(1);
        let ptr = unsafe { cuda_host_alloc(count) };
        if ptr.is_null() {
            return Err(format!("cudaHostAlloc of {count} u32 failed"));
        }
        unsafe { std::ptr::write_bytes(ptr, 0, count) };
        Ok(Self { ptr, rows, cols })
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn set(&mut self, row: usize, col: usize, value: BabyBear) {
        let i = col * self.rows + row;
        self[i] = value.value;
    }

    pub fn get(&self, row: usize, col: usize) -> BabyBear {
        BabyBear::new(self[col * self.rows + row] as u64)
    }

    /// One whole polynomial.
    pub fn column(&self, col: usize) -> &[u32] {
        &self[col * self.rows..(col + 1) * self.rows]
    }

    /// Row `row` gathered across every polynomial.
    pub fn row(&self, row: usize) -> Vec<BabyBear> {
        (0..self.cols).map(|col| self.get(row, col)).collect()
    }
}

impl Deref for PinnedSquare {
    type Target = [u32];
    fn deref(&self) -> &[u32] {
        unsafe { std::slice::from_raw_parts(self.ptr, self.rows * self.cols) }
    }
}

impl DerefMut for PinnedSquare {
    fn deref_mut(&mut self) -> &mut [u32] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr, self.rows * self.cols) }
    }
}

impl Drop for PinnedSquare {
    fn drop(&mut self) {
        unsafe { cuda_host_free(self.ptr) };
    }
}

fn run(
    values: &mut [u32],
    n: usize,
    stride: usize,
    f: unsafe extern "C" fn(*mut u32, *const u32, u32, u32, u32) -> CudaError,
    name: &str,
) -> Result<(), String> {
    assert!(n.is_power_of_two(), "NTT size must be power of 2");
    assert!(stride >= n, "stride must be at least n");
    assert!(stride > 0 && values.len() % stride == 0, "values must be whole polynomials of stride elements");
    let polys = values.len() / stride;

    // The circle: roots[i] = omega^i for i in 0..n.
    let roots: Vec<u32> = roots_of_unity(n as u32, BABYBEAR_PRIME as u32);

    let err = unsafe { f(values.as_mut_ptr(), roots.as_ptr(), n as u32, polys as u32, stride as u32) };
    if err != CUDA_SUCCESS {
        return Err(format!("{name} failed with CUDA error {err}"));
    }
    Ok(())
}

/// Forward NTT on the GPU, in place. `values` is polynomials laid end to end,
/// each `stride` long; the first `n` values of every polynomial are transformed.
/// Pass `stride == n` when the polynomials are exactly `n` long.
pub fn ntt_cuda(values: &mut [u32], n: usize, stride: usize) -> Result<(), String> {
    run(values, n, stride, cuda_ntt, "cuda_ntt")
}

/// Inverse NTT on the GPU, in place. Same layout as `ntt_cuda`.
pub fn intt_cuda(values: &mut [u32], n: usize, stride: usize) -> Result<(), String> {
    run(values, n, stride, cuda_intt, "cuda_intt")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ntt::{intt_babybear as cpu_intt, ntt_babybear as cpu_ntt};

    fn sample(n: usize) -> Vec<BabyBear> {
        (0..n).map(|i| BabyBear::new((i * 7 + 3) as u64)).collect()
    }

    fn one_poly(values: &[BabyBear]) -> PinnedSquare {
        let mut sq = PinnedSquare::new(values.len(), 1).unwrap();
        for (i, v) in values.iter().enumerate() {
            sq[i] = v.value;
        }
        sq
    }

    #[test]
    fn test_cuda_ntt_vs_cpu() {
        for log_n in [1, 4, 8, 12] {
            let n = 1usize << log_n;
            let mut cpu_values = sample(n);
            let mut gpu_values = one_poly(&cpu_values);

            cpu_ntt(&mut cpu_values);
            ntt_cuda(&mut gpu_values, n, n).unwrap();

            for (i, c) in cpu_values.iter().enumerate() {
                assert_eq!(c.value, gpu_values[i], "n={n}: mismatch at {i}");
            }
        }
    }

    #[test]
    fn test_cuda_intt_vs_cpu() {
        for log_n in [1, 4, 8, 12] {
            let n = 1usize << log_n;
            let mut cpu_values = sample(n);
            let mut gpu_values = one_poly(&cpu_values);

            cpu_intt(&mut cpu_values);
            intt_cuda(&mut gpu_values, n, n).unwrap();

            for (i, c) in cpu_values.iter().enumerate() {
                assert_eq!(c.value, gpu_values[i], "n={n}: mismatch at {i}");
            }
        }
    }

    #[test]
    fn test_cuda_roundtrip() {
        let n = 1 << 10;
        let original = sample(n);
        let mut values = one_poly(&original);
        ntt_cuda(&mut values, n, n).unwrap();
        intt_cuda(&mut values, n, n).unwrap();
        for (i, o) in original.iter().enumerate() {
            assert_eq!(o.value, values[i], "roundtrip failed at {i}");
        }
    }

    /// Many polynomials at once must equal the CPU transform of each one.
    #[test]
    fn test_cuda_many_columns_vs_cpu() {
        let (n, cols) = (1usize << 8, 37usize);
        let mut square = PinnedSquare::new(n, cols).unwrap();
        for row in 0..n {
            for col in 0..cols {
                square.set(row, col, BabyBear::new((row * 31 + col * 7 + 1) as u64));
            }
        }
        let mut cpu_columns: Vec<Vec<BabyBear>> =
            (0..cols).map(|col| (0..n).map(|row| square.get(row, col)).collect()).collect();

        ntt_cuda(&mut square, n, n).unwrap();
        for col in 0..cols {
            cpu_ntt(&mut cpu_columns[col]);
            for row in 0..n {
                assert_eq!(cpu_columns[col][row], square.get(row, col), "col={col} row={row}");
            }
        }
    }

    /// stride > n: transform only the first n values of each taller polynomial.
    #[test]
    fn test_cuda_strided_vs_cpu() {
        let (n, stride, cols) = (1usize << 6, 1usize << 8, 5usize);
        let mut square = PinnedSquare::new(stride, cols).unwrap();
        for col in 0..cols {
            for row in 0..stride {
                square.set(row, col, BabyBear::new((row * 13 + col * 3 + 2) as u64));
            }
        }
        let before: Vec<u32> = square.to_vec();

        ntt_cuda(&mut square, n, stride).unwrap();
        for col in 0..cols {
            let mut cpu: Vec<BabyBear> = (0..n).map(|row| BabyBear::new(before[col * stride + row] as u64)).collect();
            cpu_ntt(&mut cpu);
            for row in 0..stride {
                let expect = if row < n { cpu[row].value } else { before[col * stride + row] };
                assert_eq!(expect, square.get(row, col).value, "col={col} row={row}");
            }
        }
    }
}
