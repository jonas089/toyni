#include <cuda_runtime.h>
#include <stdint.h>

// GPU version of src/ntt/cpu.rs::ntt for many polynomials at once.
//
// The input is the polynomials laid end to end. Polynomial c is the n values
// starting at values[c * stride]. Usually stride == n and the array is simply
// poly 0, then poly 1, then poly 2:
//
//   values: [ p0[0] p0[1] ... p0[n-1] | p1[0] p1[1] ... p1[n-1] | p2[0] ... ]
//
// stride > n lets the caller keep each polynomial in a taller buffer and
// transform only its first n entries.

#define P          2013265921u   // the CPU's `modulus`
#define P_NEG_INV  0x77FFFFFFu   // -P^-1 mod 2^32
#define R2_MOD_P   0x45DDDDE3u   // 2^64 mod P
#define THREADS    256           // threads per block
#define CHUNK_POLYS 1024         // polynomials uploaded per chunk
#define NUM_STREAMS 3            // chunks in flight at once

// ---- field arithmetic ------------------------------------------------------
// Values live in Montgomery form (x * 2^32 mod P) so that mul needs no division.
// add/sub are the CPU's, with min() replacing the if.

__device__ inline uint32_t bb_add(uint32_t a, uint32_t b) {   // (a + b) % P
    uint32_t s = a + b;
    return min(s, s - P);
}
__device__ inline uint32_t bb_sub(uint32_t a, uint32_t b) {   // (a + P - b) % P
    uint32_t d = a - b;
    return min(d, d + P);
}
__device__ inline uint32_t mont_redc(uint32_t lo, uint32_t hi) {   // (hi:lo) / 2^32 mod P
    uint32_t m = lo * P_NEG_INV;
    uint32_t t = hi + __umulhi(m, P) + (lo != 0);
    return min(t, t - P);
}
__device__ inline uint32_t mont_mul(uint32_t a, uint32_t b) {   // (a * b) % P, Montgomery form
    return mont_redc(a * b, __umulhi(a, b));
}

// ---- elementwise kernels: one thread per value, flat over the whole chunk ----

__global__ void to_mont_kernel(uint32_t* v, uint32_t total) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < total) v[i] = mont_mul(v[i], R2_MOD_P);
}
// Leaves Montgomery form and multiplies by a plain scale (1, or n^-1 for the inverse).
__global__ void from_mont_kernel(uint32_t* v, uint32_t total, uint32_t scale) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < total) v[i] = mont_redc(v[i] * scale, __umulhi(v[i], scale));
}

// ---- per-polynomial kernels ------------------------------------------------
// Grid layout for everything below:
//   blockIdx.y                              = which polynomial
//   blockIdx.x * blockDim.x + threadIdx.x   = which position inside it
// So `poly = coeffs + blockIdx.y * n` is a plain n-element array and the rest
// of the kernel indexes it exactly like the CPU code indexes `coeffs`.

// CPU: `let mut coeffs = reverse(values);`
__global__ void bit_reverse_kernel(uint32_t* coeffs, uint32_t n, uint32_t log_n) {
    uint32_t* poly = coeffs + (size_t)blockIdx.y * n;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    uint32_t rev = __brev(i) >> (32 - log_n);
    if (i < rev) {                       // each pair is swapped once, by its lower index
        uint32_t t = poly[i]; poly[i] = poly[rev]; poly[rev] = t;
    }
}

// One stage. CPU: the body of `while len <= n`.
// A stage splits the polynomial into groups of `len` and does `half` butterflies
// per group. Butterfly b of the stage is (group b / half, position b % half):
//     even = group * len + j        CPU: lo[j]
//     odd  = even + half            CPU: hi[j]
// One thread = one butterfly of one polynomial.
__global__ void ntt_stage(uint32_t* coeffs, uint32_t n, uint32_t len, const uint32_t* roots) {
    uint32_t* poly = coeffs + (size_t)blockIdx.y * n;
    uint32_t b = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t half = len / 2;
    if (b >= n / 2) return;

    uint32_t group = b / half;
    uint32_t j     = b % half;
    uint32_t even  = group * len + j;
    uint32_t odd   = even + half;

    uint32_t w     = roots[j * (n / len)];        // CPU: let w = roots[j * stride];
    uint32_t e     = poly[even];                  // CPU: let a = lo[j];
    uint32_t w_odd = mont_mul(poly[odd], w);      // CPU: let b = hi[j] * w;
    poly[even] = bb_add(e, w_odd);                // CPU: lo[j] = a + b;
    poly[odd]  = bb_sub(e, w_odd);                // CPU: hi[j] = a - b;
}

// ---- entry points -------------------------------------------------------------
// Polynomials go through in chunks of CHUNK_POLYS; each chunk is uploaded once,
// gets one launch per stage, and is downloaded once, on its own stream. scale is
// applied on the way out. values should come from cuda_host_alloc.

static uint32_t blocks_for(uint32_t count) {   // enough blocks of THREADS to cover count
    return (count + THREADS - 1) / THREADS;
}

static cudaError_t ntt_all_polys(uint32_t* values, const uint32_t* roots, uint32_t n,
                                 uint32_t polys, uint32_t stride, uint32_t scale) {
    uint32_t log_n = 0;
    while ((1u << log_n) < n) log_n++;

    uint32_t* d_roots;
    cudaMalloc((void**)&d_roots, n * sizeof(uint32_t));
    cudaMemcpy(d_roots, roots, n * sizeof(uint32_t), cudaMemcpyHostToDevice);
    to_mont_kernel<<<blocks_for(n), THREADS>>>(d_roots, n);

    cudaStream_t streams[NUM_STREAMS];
    uint32_t* d_coeffs[NUM_STREAMS];
    for (int s = 0; s < NUM_STREAMS; s++) {
        cudaStreamCreate(&streams[s]);
        cudaMalloc((void**)&d_coeffs[s], (size_t)n * CHUNK_POLYS * sizeof(uint32_t));
    }

    size_t poly_bytes = (size_t)n * sizeof(uint32_t);
    int chunk = 0;
    for (uint32_t p0 = 0; p0 < polys; p0 += CHUNK_POLYS, chunk++) {
        uint32_t w = polys - p0 < CHUNK_POLYS ? polys - p0 : CHUNK_POLYS;   // polys in this chunk
        uint32_t total = n * w;                                             // values in this chunk
        cudaStream_t st = streams[chunk % NUM_STREAMS];
        uint32_t* coeffs = d_coeffs[chunk % NUM_STREAMS];
        uint32_t* host = values + (size_t)p0 * stride;                      // first poly of the chunk

        // Copy w polynomials of n values each. On the host consecutive polynomials
        // are `stride` values apart; on the device they are packed, n apart.
        // When stride == n this is one contiguous copy.
        cudaMemcpy2DAsync(coeffs, poly_bytes, host, (size_t)stride * sizeof(uint32_t),
                          poly_bytes, w, cudaMemcpyHostToDevice, st);

        dim3 per_value(blocks_for(n), w);        // x covers positions 0..n, y is the polynomial
        dim3 per_butterfly(blocks_for(n / 2), w); // x covers butterflies 0..n/2, y is the polynomial

        // The whole NTT for this chunk, as a sequence of kernel launches on stream
        // `st`. Launches on one stream run in order, so each step sees the previous
        // step's output; nothing here waits on the host. Every launch works on all
        // w polynomials of the chunk at once (the y dimension of the grid).

        // 1. Convert every value to Montgomery form (x -> x * 2^32 mod P) so the
        //    butterflies can use mont_mul. Flat over all n*w values.
        to_mont_kernel<<<blocks_for(total), THREADS, 0, st>>>(coeffs, total);

        // 2. Bit-reverse permutation of each polynomial's n positions.
        //    CPU: `let mut coeffs = reverse(values);`
        bit_reverse_kernel<<<per_value, THREADS, 0, st>>>(coeffs, n, log_n);

        // 3. The log2(n) butterfly stages, one launch each, len = 2, 4, ..., n.
        //    Each stage reads and writes the whole polynomial, so consecutive
        //    stages must be separate launches: the kernel boundary is the global
        //    barrier between them. CPU: the `while len <= n` loop.
        for (uint32_t len = 2; len <= n; len *= 2) {
            ntt_stage<<<per_butterfly, THREADS, 0, st>>>(coeffs, n, len, d_roots);
        }

        // 4. Leave Montgomery form and apply `scale` (1 for ntt, n^-1 for intt).
        from_mont_kernel<<<blocks_for(total), THREADS, 0, st>>>(coeffs, total, scale);

        // 5. Copy the transformed chunk back into place in `values`. Same 2D copy
        //    as the upload with source and destination swapped: device rows are n
        //    apart, host rows are `stride` apart. Queued on `st`, so it runs after
        //    step 4 and overlaps with the next chunk's work on another stream.
        cudaMemcpy2DAsync(host, (size_t)stride * sizeof(uint32_t), coeffs, poly_bytes,
                          poly_bytes, w, cudaMemcpyDeviceToHost, st);
    }

    cudaDeviceSynchronize();
    for (int s = 0; s < NUM_STREAMS; s++) {
        cudaStreamDestroy(streams[s]);
        cudaFree(d_coeffs[s]);
    }
    cudaFree(d_roots);
    return cudaGetLastError();
}

// Pinned host memory. Allocate the buffer with this once; the async copies in
// ntt_all_polys only overlap with kernels when the host buffer is pinned.
extern "C" uint32_t* cuda_host_alloc(size_t count) {
    void* p = nullptr;
    return cudaHostAlloc(&p, count * sizeof(uint32_t), cudaHostAllocDefault) == cudaSuccess ? (uint32_t*)p : nullptr;
}
extern "C" void cuda_host_free(uint32_t* p) { cudaFreeHost(p); }

extern "C" cudaError_t cuda_ntt(uint32_t* values, const uint32_t* roots, uint32_t n,
                                uint32_t polys, uint32_t stride) {
    return ntt_all_polys(values, roots, n, polys, stride, 1);
}

// CPU: intt = ntt with inv_roots, then multiply by n^-1.
extern "C" cudaError_t cuda_intt(uint32_t* values, const uint32_t* roots, uint32_t n,
                                 uint32_t polys, uint32_t stride) {
    uint32_t* inv_roots = new uint32_t[n];
    inv_roots[0] = roots[0];
    for (uint32_t i = 1; i < n; i++) inv_roots[i] = roots[n - i];

    uint64_t inv_n = 1;                       // (2^-1)^log2(n), with 2^-1 = (P+1)/2
    for (uint32_t m = n; m > 1; m /= 2) inv_n = inv_n * ((P + 1) / 2) % P;

    cudaError_t err = ntt_all_polys(values, inv_roots, n, polys, stride, (uint32_t)inv_n);
    delete[] inv_roots;
    return err;
}
