//! Metal GPU backend (`--features metal`): full prover pipeline on Apple
//! silicon.
//!
//! Every heavy stage runs as a compute kernel with results bit-identical to
//! the CPU path (verified by golden tests):
//!
//! - circle FFT butterfly layers (LDE evaluation and interpolation),
//! - zero-knowledge mask application,
//! - constraint/composition evaluation (per-AIR MSL generated at runtime
//!   from [`crate::air::Air::metal_transitions`]),
//! - boundary quotient terms,
//! - DEEP quotient batching,
//! - FRI folds,
//! - SHA-256 Merkle tree construction.
//!
//! Buffers use shared storage (Apple unified memory), so "upload" and
//! "readback" are cheap. The context is cached per thread; set
//! `CIRCLE_STARK_BACKEND=cpu` to force the CPU path in a metal-enabled build.

use std::cell::RefCell;
use std::collections::HashMap;
use std::ffi::c_void;

use metal::{
    Buffer, CommandQueue, CompileOptions, ComputePipelineState, Device, MTLResourceOptions,
    MTLSize,
};

use crate::engine::circle::cfft::Twiddles;
use crate::field::{M31, QM31};
use crate::merkle::{Hash, MerkleTree};

const KERNELS_SRC: &str = include_str!("kernels.metal");

/// Minimum problem size before dispatching to the GPU (below this, kernel
/// launch overhead dominates).
pub const GPU_THRESHOLD: usize = 1 << 14;

/// A boundary constraint prepared for the GPU composition pass.
#[derive(Debug, Clone, Copy)]
pub struct GpuBoundary {
    pub px: M31,
    pub py: M31,
    pub value: M31,
    pub col: usize,
    pub beta_pow: QM31,
}

pub struct MetalContext {
    device: Device,
    queue: CommandQueue,
    pipelines: HashMap<&'static str, ComputePipelineState>,
    /// Runtime-compiled composition pipelines, keyed by AIR label.
    air_pipelines: HashMap<String, ComputePipelineState>,
}

thread_local! {
    static CONTEXT: RefCell<Option<Option<MetalContext>>> = const { RefCell::new(None) };
}

/// Run `f` with the thread's Metal context, initializing it on first use.
/// Returns `None` if Metal is unavailable or disabled via env.
pub fn with_context<R>(f: impl FnOnce(&mut MetalContext) -> R) -> Option<R> {
    if std::env::var("CIRCLE_STARK_BACKEND").as_deref() == Ok("cpu") {
        return None;
    }
    CONTEXT.with(|slot| {
        let mut slot = slot.borrow_mut();
        if slot.is_none() {
            *slot = Some(MetalContext::new());
        }
        slot.as_mut().unwrap().as_mut().map(f)
    })
}

/// True if the GPU path is available on this thread.
pub fn available() -> bool {
    with_context(|_| ()).is_some()
}

const BASE_KERNELS: &[&str] = &[
    "cfft_fwd_layer_m31",
    "cfft_inv_layer_m31",
    "m31_scale",
    "mask_add",
    "boundary_term",
    "deep_combine",
    "fri_fold",
    "merkle_leaves",
    "merkle_nodes",
];

impl MetalContext {
    fn new() -> Option<Self> {
        let device = Device::system_default()?;
        let queue = device.new_command_queue();
        let library = device
            .new_library_with_source(KERNELS_SRC, &CompileOptions::new())
            .map_err(|e| eprintln!("circle-stark: metal shader compile failed: {e}"))
            .ok()?;
        let mut pipelines = HashMap::new();
        for name in BASE_KERNELS {
            let f = library.get_function(name, None).ok()?;
            let pso = device.new_compute_pipeline_state_with_function(&f).ok()?;
            pipelines.insert(*name, pso);
        }
        Some(Self { device, queue, pipelines, air_pipelines: HashMap::new() })
    }

    // ── buffer helpers ─────────────────────────────────────────────────

    fn buffer_from<T: Copy>(&self, data: &[T]) -> Buffer {
        self.device.new_buffer_with_data(
            data.as_ptr() as *const c_void,
            std::mem::size_of_val(data) as u64,
            MTLResourceOptions::StorageModeShared,
        )
    }

    fn empty_buffer(&self, bytes: usize) -> Buffer {
        self.device
            .new_buffer(bytes as u64, MTLResourceOptions::StorageModeShared)
    }

    fn read_buffer<T: Copy>(buf: &Buffer, len: usize) -> Vec<T> {
        unsafe { std::slice::from_raw_parts(buf.contents() as *const T, len).to_vec() }
    }

    /// Encode one pipeline dispatch over `n` threads.
    fn dispatch(
        &self,
        encoder: &metal::ComputeCommandEncoderRef,
        pso: &ComputePipelineState,
        buffers: &[&Buffer],
        bytes: Option<(&[u8], u64)>,
        n: usize,
    ) {
        encoder.set_compute_pipeline_state(pso);
        for (i, b) in buffers.iter().enumerate() {
            encoder.set_buffer(i as u64, Some(b), 0);
        }
        if let Some((data, index)) = bytes {
            encoder.set_bytes(index, data.len() as u64, data.as_ptr() as *const c_void);
        }
        let tg = 256.min(n) as u64;
        encoder.dispatch_threads(MTLSize::new(n as u64, 1, 1), MTLSize::new(tg, 1, 1));
    }

    fn run(&self, encode: impl FnOnce(&metal::ComputeCommandEncoderRef)) {
        objc::rc::autoreleasepool(|| {
            let cb = self.queue.new_command_buffer();
            let encoder = cb.new_compute_command_encoder();
            encode(encoder);
            encoder.end_encoding();
            cb.commit();
            cb.wait_until_completed();
        })
    }

    // ── CFFT ───────────────────────────────────────────────────────────

    /// In-place forward CFFT (coefficients → values) on the GPU.
    pub fn evaluate_m31(&self, coeffs: &mut [M31], twiddles: &Twiddles) {
        let n = twiddles.log_size as usize;
        assert_eq!(coeffs.len(), 1 << n);
        let size = coeffs.len();
        let data = self.buffer_from(coeffs);
        let tw_buffers: Vec<Buffer> =
            twiddles.forward.iter().map(|l| self.buffer_from(l)).collect();

        self.run(|enc| {
            let pso = &self.pipelines["cfft_fwd_layer_m31"];
            let mut step: u32 = 1;
            for l in (1..n).rev() {
                self.dispatch(
                    enc,
                    pso,
                    &[&data, &tw_buffers[l]],
                    Some((&step.to_le_bytes(), 2)),
                    size / 2,
                );
                step *= 2;
            }
            let half_step = (size as u32) / 2;
            self.dispatch(
                enc,
                pso,
                &[&data, &tw_buffers[0]],
                Some((&half_step.to_le_bytes(), 2)),
                size / 2,
            );
        });
        coeffs.copy_from_slice(&Self::read_buffer(&data, size));
    }

    /// In-place inverse CFFT (values → coefficients) on the GPU.
    pub fn interpolate_m31(&self, values: &mut [M31], twiddles: &Twiddles) {
        let n = twiddles.log_size as usize;
        assert_eq!(values.len(), 1 << n);
        let size = values.len();
        let data = self.buffer_from(values);
        let tw_buffers: Vec<Buffer> =
            twiddles.inverse.iter().map(|l| self.buffer_from(l)).collect();
        let scale = M31::TWO.pow(n as u64).inverse();

        self.run(|enc| {
            let pso = &self.pipelines["cfft_inv_layer_m31"];
            let half_step = (size as u32) / 2;
            self.dispatch(
                enc,
                pso,
                &[&data, &tw_buffers[0]],
                Some((&half_step.to_le_bytes(), 2)),
                size / 2,
            );
            let mut step: u32 = (size as u32) / 4;
            for tw in tw_buffers.iter().take(n).skip(1) {
                self.dispatch(
                    enc,
                    pso,
                    &[&data, tw],
                    Some((&step.to_le_bytes(), 2)),
                    size / 2,
                );
                step /= 2;
            }
            self.dispatch(
                enc,
                &self.pipelines["m31_scale"],
                &[&data],
                Some((&scale.0.to_le_bytes(), 1)),
                size,
            );
        });
        values.copy_from_slice(&Self::read_buffer(&data, size));
    }

    /// QM31 transforms decompose into four base-coordinate M31 transforms
    /// (the twiddles are base-field, so the CFFT is M31-linear).
    pub fn interpolate_qm31(&self, values: &mut [QM31], twiddles: &Twiddles) {
        let len = values.len();
        let mut coords: [Vec<M31>; 4] = std::array::from_fn(|_| vec![M31::ZERO; len]);
        for (i, v) in values.iter().enumerate() {
            let arr = v.to_m31_array();
            for k in 0..4 {
                coords[k][i] = arr[k];
            }
        }
        for coord in coords.iter_mut() {
            self.interpolate_m31(coord, twiddles);
        }
        for (i, v) in values.iter_mut().enumerate() {
            *v = QM31::from_m31_array([coords[0][i], coords[1][i], coords[2][i], coords[3][i]]);
        }
    }

    /// lde[i] += vn[i] · r[i]  (zero-knowledge masking).
    pub fn mask_add(&self, lde: &mut [M31], vn: &[M31], r: &[M31]) {
        let n = lde.len();
        let data = self.buffer_from(lde);
        let vn_buf = self.buffer_from(vn);
        let r_buf = self.buffer_from(r);
        self.run(|enc| {
            self.dispatch(enc, &self.pipelines["mask_add"], &[&data, &vn_buf, &r_buf], None, n);
        });
        lde.copy_from_slice(&Self::read_buffer(&data, n));
    }

    // ── FRI fold ───────────────────────────────────────────────────────

    pub fn fri_fold(&self, values: &[QM31], inv_twiddles: &[M31], lambda: QM31) -> Vec<QM31> {
        let half = values.len() / 2;
        assert_eq!(inv_twiddles.len(), half);
        let input = self.buffer_from(values);
        let tw = self.buffer_from(inv_twiddles);
        let out = self.empty_buffer(half * 16);

        #[repr(C)]
        #[derive(Clone, Copy)]
        struct FoldParams {
            half_size: u32,
            _pad: [u32; 3],
            lambda: [u32; 4],
        }
        let params = FoldParams {
            half_size: half as u32,
            _pad: [0; 3],
            lambda: lambda.to_m31_array().map(|m| m.0),
        };
        let bytes = unsafe {
            std::slice::from_raw_parts(
                &params as *const FoldParams as *const u8,
                std::mem::size_of::<FoldParams>(),
            )
        };
        self.run(|enc| {
            self.dispatch(enc, &self.pipelines["fri_fold"], &[&out, &input, &tw], Some((bytes, 3)), half);
        });
        Self::read_buffer(&out, half)
    }

    // ── Merkle ─────────────────────────────────────────────────────────

    /// Build a Merkle tree over fixed-size leaves on the GPU.
    pub fn merkle(&self, packed_leaves: &[u8], leaf_len: usize, count: usize) -> MerkleTree {
        assert!(count.is_power_of_two());
        assert!(leaf_len < 159, "leaf too large for the GPU hasher");
        assert_eq!(packed_leaves.len(), leaf_len * count);
        let data = self.buffer_from(packed_leaves);

        // Level 0: leaf hashes.
        let mut level_buffers: Vec<Buffer> = vec![self.empty_buffer(count * 32)];
        let mut sizes = vec![count];
        while *sizes.last().unwrap() > 1 {
            let next = sizes.last().unwrap() / 2;
            level_buffers.push(self.empty_buffer(next * 32));
            sizes.push(next);
        }

        #[repr(C)]
        struct HashParams {
            stride: u32,
            len: u32,
        }
        let params = HashParams { stride: leaf_len as u32, len: leaf_len as u32 };
        let pbytes = unsafe {
            std::slice::from_raw_parts(&params as *const HashParams as *const u8, 8)
        };

        self.run(|enc| {
            self.dispatch(
                enc,
                &self.pipelines["merkle_leaves"],
                &[&level_buffers[0], &data],
                Some((pbytes, 2)),
                count,
            );
            for (l, &level_size) in sizes.iter().enumerate().skip(1) {
                let (prev, cur) = level_buffers.split_at(l);
                self.dispatch(
                    enc,
                    &self.pipelines["merkle_nodes"],
                    &[&cur[0], prev.last().unwrap()],
                    None,
                    level_size,
                );
            }
        });

        let levels: Vec<Vec<Hash>> = level_buffers
            .iter()
            .zip(&sizes)
            .map(|(buf, &s)| Self::read_buffer::<Hash>(buf, s))
            .collect();
        MerkleTree { levels }
    }

    // ── Composition (AIR-templated) ────────────────────────────────────

    /// Get or compile the composition pipeline for an AIR.
    fn air_pipeline(
        &mut self,
        label: &str,
        air_msl: &str,
        num_offsets: usize,
        num_cols: usize,
        num_constraints: usize,
    ) -> Option<ComputePipelineState> {
        if let Some(pso) = self.air_pipelines.get(label) {
            return Some(pso.clone());
        }
        let src = format!(
            "{KERNELS_SRC}\n\
             #define NUM_OFFSETS {num_offsets}\n\
             #define NUM_COLS {num_cols}\n\
             #define NUM_CONSTRAINTS {num_constraints}\n\
             {air_msl}\n\
             {COMPOSITION_TEMPLATE}"
        );
        let library = self
            .device
            .new_library_with_source(&src, &CompileOptions::new())
            .map_err(|e| eprintln!("circle-stark: AIR shader compile failed for {label}: {e}"))
            .ok()?;
        let f = library.get_function("composition_transitions", None).ok()?;
        let pso = self.device.new_compute_pipeline_state_with_function(&f).ok()?;
        self.air_pipelines.insert(label.to_string(), pso.clone());
        Some(pso)
    }

    /// The full composition polynomial (transitions + boundary quotients)
    /// over the whole domain, in one command buffer with one readback.
    #[allow(clippy::too_many_arguments)]
    pub fn composition(
        &mut self,
        label: &str,
        air_msl: &str,
        trace_ldes: &[Vec<M31>],
        inv_vn: &[M31],
        selector: Option<&[M31]>,
        beta_pows: &[QM31],
        num_offsets: usize,
        rotation_step: usize,
        points_x: &[M31],
        points_y: &[M31],
        boundaries: &[GpuBoundary],
    ) -> Option<Vec<QM31>> {
        let num_cols = trace_ldes.len();
        let domain_size = trace_ldes[0].len();
        let pso = self.air_pipeline(label, air_msl, num_offsets, num_cols, beta_pows.len())?;

        // Column-major packed trace, shared by both kernels.
        let mut packed = Vec::with_capacity(num_cols * domain_size);
        for col in trace_ldes {
            packed.extend_from_slice(col);
        }
        let trace_buf = self.buffer_from(&packed);
        let vn_buf = self.buffer_from(inv_vn);
        let one = vec![M31::ONE; 1];
        let sel_buf = self.buffer_from(selector.unwrap_or(&one));
        let beta_flat: Vec<[u32; 4]> =
            beta_pows.iter().map(|b| b.to_m31_array().map(|m| m.0)).collect();
        let beta_buf = self.buffer_from(&beta_flat);
        let x_buf = self.buffer_from(points_x);
        let y_buf = self.buffer_from(points_y);
        let out = self.empty_buffer(domain_size * 16);

        #[repr(C)]
        struct CompParams {
            domain_size: u32,
            half: u32,
            rot_step: u32,
            use_selector: u32,
        }
        let params = CompParams {
            domain_size: domain_size as u32,
            half: (domain_size / 2) as u32,
            rot_step: rotation_step as u32,
            use_selector: selector.is_some() as u32,
        };
        let pbytes = unsafe {
            std::slice::from_raw_parts(&params as *const CompParams as *const u8, 16)
        };

        #[repr(C)]
        struct BoundaryParams {
            px: u32,
            py: u32,
            value: u32,
            col: u32,
            beta_pow: [u32; 4],
            domain_size: u32,
            _pad: [u32; 3],
        }
        let boundary_params: Vec<BoundaryParams> = boundaries
            .iter()
            .map(|b| BoundaryParams {
                px: b.px.0,
                py: b.py.0,
                value: b.value.0,
                col: b.col as u32,
                beta_pow: b.beta_pow.to_m31_array().map(|m| m.0),
                domain_size: domain_size as u32,
                _pad: [0; 3],
            })
            .collect();

        self.run(|enc| {
            self.dispatch(
                enc,
                &pso,
                &[&out, &trace_buf, &vn_buf, &sel_buf, &beta_buf],
                Some((pbytes, 5)),
                domain_size,
            );
            for bp in &boundary_params {
                let bbytes = unsafe {
                    std::slice::from_raw_parts(
                        bp as *const BoundaryParams as *const u8,
                        std::mem::size_of::<BoundaryParams>(),
                    )
                };
                self.dispatch(
                    enc,
                    &self.pipelines["boundary_term"],
                    &[&out, &trace_buf, &x_buf, &y_buf],
                    Some((bbytes, 4)),
                    domain_size,
                );
            }
        });
        Some(Self::read_buffer(&out, domain_size))
    }

    /// The DEEP batched word over the whole domain.
    #[allow(clippy::too_many_arguments)]
    pub fn deep_combine(
        &self,
        trace_ldes: &[Vec<M31>],
        composition: &[QM31],
        points_x: &[M31],
        points_y: &[M31],
        mask_points: &[crate::engine::circle::geometry::CirclePoint<QM31>],
        ood_trace: &[Vec<QM31>],
        ood_composition: QM31,
        mu_pows: &[QM31],
    ) -> Vec<QM31> {
        let domain_size = composition.len();
        let num_cols = trace_ldes.len();
        let num_offsets = mask_points.len();

        let mut packed = Vec::with_capacity(num_cols * domain_size);
        for c in trace_ldes {
            packed.extend_from_slice(c);
        }
        let trace_buf = self.buffer_from(&packed);
        let comp_buf = self.buffer_from(composition);
        let x_buf = self.buffer_from(points_x);
        let y_buf = self.buffer_from(points_y);

        let mut mask_flat: Vec<[u32; 4]> = Vec::with_capacity(num_offsets * 2);
        for z in mask_points {
            mask_flat.push(z.x.to_m31_array().map(|m| m.0));
            mask_flat.push(z.y.to_m31_array().map(|m| m.0));
        }
        let mask_buf = self.buffer_from(&mask_flat);

        let ood_flat: Vec<[u32; 4]> = ood_trace
            .iter()
            .flat_map(|row| row.iter().map(|v| v.to_m31_array().map(|m| m.0)))
            .collect();
        let ood_buf = self.buffer_from(&ood_flat);
        let mu_flat: Vec<[u32; 4]> = mu_pows.iter().map(|v| v.to_m31_array().map(|m| m.0)).collect();
        let mu_buf = self.buffer_from(&mu_flat);
        let out = self.empty_buffer(domain_size * 16);

        #[repr(C)]
        struct DeepParams {
            domain_size: u32,
            num_cols: u32,
            num_offsets: u32,
            _pad: u32,
            ood_comp: [u32; 4],
        }
        let params = DeepParams {
            domain_size: domain_size as u32,
            num_cols: num_cols as u32,
            num_offsets: num_offsets as u32,
            _pad: 0,
            ood_comp: ood_composition.to_m31_array().map(|m| m.0),
        };
        let pbytes = unsafe {
            std::slice::from_raw_parts(&params as *const DeepParams as *const u8, 32)
        };
        self.run(|enc| {
            self.dispatch(
                enc,
                &self.pipelines["deep_combine"],
                &[&out, &trace_buf, &comp_buf, &x_buf, &y_buf, &mask_buf, &ood_buf, &mu_buf],
                Some((pbytes, 8)),
                domain_size,
            );
        });
        Self::read_buffer(&out, domain_size)
    }
}

/// The AIR-agnostic composition kernel; `air_transitions` is supplied by the
/// AIR (see [`crate::air::Air::metal_transitions`]).
const COMPOSITION_TEMPLATE: &str = r#"
struct CompParams {
    uint domain_size;
    uint half_size;
    uint rot_step;
    uint use_selector;
};

inline uint rot_idx(uint i, uint steps, uint half_size) {
    if (i < half_size) { return (i + steps) % half_size; }
    uint j = i - half_size;
    return half_size + (j + half_size - (steps % half_size)) % half_size;
}

kernel void composition_transitions(
    device uint4* out [[buffer(0)]],
    const device uint* trace [[buffer(1)]],
    const device uint* inv_vn [[buffer(2)]],
    const device uint* selector [[buffer(3)]],
    const device uint4* beta_pows [[buffer(4)]],
    constant CompParams& p [[buffer(5)]],
    uint i [[thread_position_in_grid]])
{
    uint m[NUM_OFFSETS][NUM_COLS];
    for (uint k = 0; k < NUM_OFFSETS; k++) {
        uint idx = rot_idx(i, k * p.rot_step, p.half_size);
        for (uint c = 0; c < NUM_COLS; c++) {
            m[k][c] = trace[c * p.domain_size + idx];
        }
    }
    uint cst[NUM_CONSTRAINTS];
    air_transitions(m, cst);
    uint scale = m31_mul(p.use_selector != 0 ? selector[i] : 1u, inv_vn[i]);
    uint4 acc = uint4(0);
    for (uint t = 0; t < NUM_CONSTRAINTS; t++) {
        acc = qm31_add(acc, qm31_mul_m31(beta_pows[t], m31_mul(cst[t], scale)));
    }
    out[i] = acc;
}
"#;
