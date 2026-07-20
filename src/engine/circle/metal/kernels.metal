// Metal compute kernels for the circle STARK prover (M3-class Apple GPUs).
//
// Field layout conventions (little-endian, matching the Rust structs):
//   M31   = uint            (reduced, < P)
//   CM31  = uint2  (a, b)   : a + b·i
//   QM31  = uint4  (aa, ab, ba, bb) : (aa + ab·i) + (ba + bb·i)·u
//
// All kernels are pure functions of their buffers; results are bit-identical
// to the CPU implementation (verified by golden tests).

#include <metal_stdlib>
using namespace metal;

constant uint P = 0x7fffffffu;

// ── M31 ────────────────────────────────────────────────────────────────────

inline uint m31_reduce(ulong x) {
    ulong first = (x >> 31) + (x & P);
    ulong second = (first >> 31) + (first & P);
    uint r = (uint)second;
    return r >= P ? r - P : r;
}

inline uint m31_add(uint a, uint b) {
    uint s = a + b;
    return s >= P ? s - P : s;
}

inline uint m31_sub(uint a, uint b) {
    return a >= b ? a - b : a + (P - b);
}

inline uint m31_neg(uint a) {
    return a == 0 ? 0 : P - a;
}

inline uint m31_mul(uint a, uint b) {
    return m31_reduce((ulong)a * (ulong)b);
}

inline uint m31_pow(uint a, uint e) {
    uint r = 1;
    while (e != 0) {
        if (e & 1) { r = m31_mul(r, a); }
        a = m31_mul(a, a);
        e >>= 1;
    }
    return r;
}

inline uint m31_inv(uint a) {
    return m31_pow(a, P - 2);
}

// ── CM31 ───────────────────────────────────────────────────────────────────

inline uint2 cm31_add(uint2 a, uint2 b) { return uint2(m31_add(a.x, b.x), m31_add(a.y, b.y)); }
inline uint2 cm31_sub(uint2 a, uint2 b) { return uint2(m31_sub(a.x, b.x), m31_sub(a.y, b.y)); }

inline uint2 cm31_mul(uint2 a, uint2 b) {
    return uint2(
        m31_sub(m31_mul(a.x, b.x), m31_mul(a.y, b.y)),
        m31_add(m31_mul(a.x, b.y), m31_mul(a.y, b.x)));
}

inline uint2 cm31_mul_m31(uint2 a, uint s) { return uint2(m31_mul(a.x, s), m31_mul(a.y, s)); }

inline uint2 cm31_inv(uint2 a) {
    uint norm = m31_add(m31_mul(a.x, a.x), m31_mul(a.y, a.y));
    uint ninv = m31_inv(norm);
    return uint2(m31_mul(a.x, ninv), m31_neg(m31_mul(a.y, ninv)));
}

// ── QM31 (u^2 = 2 + i) ─────────────────────────────────────────────────────

constant uint2 QR = uint2(2u, 1u);

inline uint4 qm31_add(uint4 a, uint4 b) {
    return uint4(cm31_add(a.xy, b.xy), cm31_add(a.zw, b.zw));
}

inline uint4 qm31_sub(uint4 a, uint4 b) {
    return uint4(cm31_sub(a.xy, b.xy), cm31_sub(a.zw, b.zw));
}

inline uint4 qm31_mul(uint4 a, uint4 b) {
    uint2 lo = cm31_add(cm31_mul(a.xy, b.xy), cm31_mul(QR, cm31_mul(a.zw, b.zw)));
    uint2 hi = cm31_add(cm31_mul(a.xy, b.zw), cm31_mul(a.zw, b.xy));
    return uint4(lo, hi);
}

inline uint4 qm31_mul_m31(uint4 a, uint s) {
    return uint4(cm31_mul_m31(a.xy, s), cm31_mul_m31(a.zw, s));
}

inline uint4 qm31_mul_cm31(uint4 a, uint2 s) {
    return uint4(cm31_mul(a.xy, s), cm31_mul(a.zw, s));
}

inline uint4 qm31_inv(uint4 a) {
    // (x + yu)^-1 = (x - yu) / (x^2 - R y^2)
    uint2 denom = cm31_sub(cm31_mul(a.xy, a.xy), cm31_mul(QR, cm31_mul(a.zw, a.zw)));
    uint2 dinv = cm31_inv(denom);
    return uint4(cm31_mul(a.xy, dinv), cm31_mul(uint2(m31_neg(a.z), m31_neg(a.w)), dinv));
}

// ── CFFT butterfly layers ──────────────────────────────────────────────────
// One dispatch per layer; thread t handles pair (i, i + step) where
// i = (t / step)·2·step + (t % step), twiddle index t % step.

struct LayerParams {
    uint step;
};

kernel void cfft_fwd_layer_m31(
    device uint* data [[buffer(0)]],
    const device uint* twiddles [[buffer(1)]],
    constant LayerParams& p [[buffer(2)]],
    uint t [[thread_position_in_grid]])
{
    uint s = p.step;
    uint k = t % s;
    uint i = (t / s) * 2 * s + k;
    uint a = data[i];
    uint tb = m31_mul(twiddles[k], data[i + s]);
    data[i] = m31_add(a, tb);
    data[i + s] = m31_sub(a, tb);
}

kernel void cfft_inv_layer_m31(
    device uint* data [[buffer(0)]],
    const device uint* inv_twiddles [[buffer(1)]],
    constant LayerParams& p [[buffer(2)]],
    uint t [[thread_position_in_grid]])
{
    uint s = p.step;
    uint k = t % s;
    uint i = (t / s) * 2 * s + k;
    uint a = data[i];
    uint b = data[i + s];
    data[i] = m31_add(a, b);
    data[i + s] = m31_mul(m31_sub(a, b), inv_twiddles[k]);
}

kernel void m31_scale(
    device uint* data [[buffer(0)]],
    constant uint& s [[buffer(1)]],
    uint i [[thread_position_in_grid]])
{
    data[i] = m31_mul(data[i], s);
}

// Add the zero-knowledge mask: lde[i] += vn[i] * r_lde[i].
kernel void mask_add(
    device uint* lde [[buffer(0)]],
    const device uint* vn [[buffer(1)]],
    const device uint* r_lde [[buffer(2)]],
    uint i [[thread_position_in_grid]])
{
    lde[i] = m31_add(lde[i], m31_mul(vn[i], r_lde[i]));
}

// ── Composition: boundary quotient terms ───────────────────────────────────
// acc[i] += beta_pow · (trace_col[i] - value) / v_Pb(P_i), with
// v_Pb(P) = 1 - (P·Pb^{-1}).x - i·(P·Pb^{-1}).y  ∈ CM31.

struct BoundaryParams {
    uint px;
    uint py;
    uint value;
    uint col;
    uint4 beta_pow;
    uint domain_size;
    uint _pad0;
    uint _pad1;
    uint _pad2;
};

kernel void boundary_term(
    device uint4* acc [[buffer(0)]],
    const device uint* trace [[buffer(1)]],       // column-major [col][i]
    const device uint* points_x [[buffer(2)]],
    const device uint* points_y [[buffer(3)]],
    constant BoundaryParams& b [[buffer(4)]],
    uint i [[thread_position_in_grid]])
{
    uint x = points_x[i];
    uint y = points_y[i];
    // rel = P · J(Pb): (x·bx + y·by, y·bx - x·by)
    uint rx = m31_add(m31_mul(x, b.px), m31_mul(y, b.py));
    uint ry = m31_sub(m31_mul(y, b.px), m31_mul(x, b.py));
    uint2 v = uint2(m31_sub(1u, rx), m31_neg(ry));
    uint2 vinv = cm31_inv(v);
    uint diff = m31_sub(trace[b.col * b.domain_size + i], b.value);
    uint2 q = cm31_mul_m31(vinv, diff);
    acc[i] = qm31_add(acc[i], qm31_mul_cm31(b.beta_pow, q));
}

// ── DEEP quotient batching ─────────────────────────────────────────────────
// u[i] = Σ_{k,c} µ^(kW+c)·(trace_c[i] - ood[k][c])·inv(v_{z_k}(P_i))
//      + µ^(KW)·(comp[i] - ood_comp)·inv(v_{z_0}(P_i))

struct DeepParams {
    uint domain_size;
    uint num_cols;
    uint num_offsets;
    uint _pad;
    uint4 ood_comp;
};

kernel void deep_combine(
    device uint4* out [[buffer(0)]],
    const device uint* trace [[buffer(1)]],       // column-major
    const device uint4* comp [[buffer(2)]],
    const device uint* points_x [[buffer(3)]],
    const device uint* points_y [[buffer(4)]],
    const device uint4* mask_pts [[buffer(5)]],   // per offset: [x, y] pairs (2 uint4 each)
    const device uint4* ood [[buffer(6)]],        // [k*num_cols + c]
    const device uint4* mu_pows [[buffer(7)]],    // num_offsets*num_cols + 1
    constant DeepParams& p [[buffer(8)]],
    uint i [[thread_position_in_grid]])
{
    uint x = points_x[i];
    uint y = points_y[i];
    uint4 acc = uint4(0);
    uint term = 0;
    for (uint k = 0; k < p.num_offsets; k++) {
        uint4 zx = mask_pts[2 * k];
        uint4 zy = mask_pts[2 * k + 1];
        // rel = P · J(z): x-coord = x·zx + y·zy ; y-coord = y·zx - x·zy
        uint4 rx = qm31_add(qm31_mul_m31(zx, x), qm31_mul_m31(zy, y));
        uint4 ry = qm31_sub(qm31_mul_m31(zx, y), qm31_mul_m31(zy, x));
        // v = 1 - rx - i·ry
        uint4 v = uint4(m31_sub(1u, rx.x), rx.y, rx.z, rx.w);
        v.y = m31_neg(v.y); v.z = m31_neg(v.z); v.w = m31_neg(v.w);
        // subtract i·ry: i·(a+bi) = -b + ai (per CM31 coordinate pair)
        uint4 iry = uint4(m31_neg(ry.y), ry.x, m31_neg(ry.w), ry.z);
        v = qm31_sub(v, iry);
        uint4 vinv = qm31_inv(v);
        if (k == 0) {
            uint4 cdiff = qm31_sub(comp[i], p.ood_comp);
            uint compterm = p.num_offsets * p.num_cols;
            out[i] = qm31_mul(qm31_mul(mu_pows[compterm], cdiff), vinv);
        }
        for (uint c = 0; c < p.num_cols; c++) {
            uint tv = trace[c * p.domain_size + i];
            uint4 diff = qm31_sub(uint4(tv, 0, 0, 0), ood[k * p.num_cols + c]);
            acc = qm31_add(acc, qm31_mul(qm31_mul(mu_pows[term], diff), vinv));
            term += 1;
        }
    }
    out[i] = qm31_add(out[i], acc);
}

// ── FRI ────────────────────────────────────────────────────────────────────

struct FoldParams {
    uint half_size;
    uint _pad0;
    uint _pad1;
    uint _pad2;
    uint4 lambda;
};

// out[i] = (a+b)/2 + λ·(a-b)/(2t_i),  a = in[i], b = in[i+half]
kernel void fri_fold(
    device uint4* out [[buffer(0)]],
    const device uint4* in [[buffer(1)]],
    const device uint* inv_twiddles [[buffer(2)]],
    constant FoldParams& p [[buffer(3)]],
    uint i [[thread_position_in_grid]])
{
    const uint HALF = 0x40000000u; // 1/2 mod P
    uint4 a = in[i];
    uint4 b = in[i + p.half_size];
    uint4 avg = qm31_mul_m31(qm31_add(a, b), HALF);
    uint4 diff = qm31_mul_m31(qm31_sub(a, b), m31_mul(HALF, inv_twiddles[i]));
    out[i] = qm31_add(avg, qm31_mul(p.lambda, diff));
}

// ── SHA-256 (Merkle hashing) ───────────────────────────────────────────────

constant uint SHA_K[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
    0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
    0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
    0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
    0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
    0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

struct Sha256State {
    uint h[8];
};

inline uint rotr(uint x, uint n) { return (x >> n) | (x << (32 - n)); }

inline void sha256_init(thread Sha256State& s) {
    s.h[0] = 0x6a09e667u; s.h[1] = 0xbb67ae85u; s.h[2] = 0x3c6ef372u; s.h[3] = 0xa54ff53au;
    s.h[4] = 0x510e527fu; s.h[5] = 0x9b05688cu; s.h[6] = 0x1f83d9abu; s.h[7] = 0x5be0cd19u;
}

inline void sha256_block(thread Sha256State& s, thread const uchar* block) {
    uint w[64];
    for (uint i = 0; i < 16; i++) {
        w[i] = ((uint)block[4 * i] << 24) | ((uint)block[4 * i + 1] << 16)
             | ((uint)block[4 * i + 2] << 8) | (uint)block[4 * i + 3];
    }
    for (uint i = 16; i < 64; i++) {
        uint s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
        uint s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    uint a = s.h[0], b = s.h[1], c = s.h[2], d = s.h[3];
    uint e = s.h[4], f = s.h[5], g = s.h[6], h = s.h[7];
    for (uint i = 0; i < 64; i++) {
        uint S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        uint ch = (e & f) ^ (~e & g);
        uint t1 = h + S1 + ch + SHA_K[i] + w[i];
        uint S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        uint maj = (a & b) ^ (a & c) ^ (b & c);
        uint t2 = S0 + maj;
        h = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    s.h[0] += a; s.h[1] += b; s.h[2] += c; s.h[3] += d;
    s.h[4] += e; s.h[5] += f; s.h[6] += g; s.h[7] += h;
}

// Hash `msg_len` bytes (already tagged) into a 32-byte digest.
// MAX_MSG must cover tag + leaf payload; leaves here are ≤ 128 bytes.
constant uint MAX_MSG = 160;

inline void sha256_digest(thread const uchar* msg, uint msg_len, device uchar* out) {
    Sha256State s;
    sha256_init(s);
    // Process full blocks.
    uint off = 0;
    while (off + 64 <= msg_len) {
        sha256_block(s, msg + off);
        off += 64;
    }
    // Final block(s) with padding.
    uchar tail[128];
    uint rem = msg_len - off;
    for (uint i = 0; i < rem; i++) { tail[i] = msg[off + i]; }
    tail[rem] = 0x80;
    uint tail_len = (rem + 9 <= 64) ? 64 : 128;
    for (uint i = rem + 1; i < tail_len - 8; i++) { tail[i] = 0; }
    ulong bitlen = (ulong)msg_len * 8;
    for (uint i = 0; i < 8; i++) {
        tail[tail_len - 1 - i] = (uchar)(bitlen >> (8 * i));
    }
    sha256_block(s, tail);
    if (tail_len == 128) { sha256_block(s, tail + 64); }
    for (uint i = 0; i < 8; i++) {
        out[4 * i] = (uchar)(s.h[i] >> 24);
        out[4 * i + 1] = (uchar)(s.h[i] >> 16);
        out[4 * i + 2] = (uchar)(s.h[i] >> 8);
        out[4 * i + 3] = (uchar)(s.h[i]);
    }
}

struct HashParams {
    uint stride;   // bytes between consecutive leaves in the input buffer
    uint len;      // payload bytes per leaf
};

// Leaf hash: out[i] = SHA256(0x00 || data[i·stride .. i·stride+len])
kernel void merkle_leaves(
    device uchar* out [[buffer(0)]],
    const device uchar* data [[buffer(1)]],
    constant HashParams& p [[buffer(2)]],
    uint i [[thread_position_in_grid]])
{
    uchar msg[MAX_MSG];
    msg[0] = 0x00;
    for (uint j = 0; j < p.len; j++) {
        msg[1 + j] = data[i * p.stride + j];
    }
    sha256_digest(msg, p.len + 1, out + 32 * i);
}

// Node hash: out[i] = SHA256(0x01 || in[2i] || in[2i+1])
kernel void merkle_nodes(
    device uchar* out [[buffer(0)]],
    const device uchar* in [[buffer(1)]],
    uint i [[thread_position_in_grid]])
{
    uchar msg[65];
    msg[0] = 0x01;
    for (uint j = 0; j < 64; j++) {
        msg[1 + j] = in[64 * i + j];
    }
    sha256_digest(msg, 65, out + 32 * i);
}
