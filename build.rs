use std::env;
use std::path::PathBuf;
use std::process::Command;

fn has_nvcc() -> bool {
    Command::new("nvcc")
        .arg("--version")
        .output()
        .is_ok()
}

/// Returns true if the installed nvcc lists `arch` (e.g. "compute_120") in
/// its supported arch list. Used to gate Blackwell (RTX 5090) codegen on
/// toolkits old enough to predate sm_120.
fn nvcc_supports_arch(arch: &str) -> bool {
    let out = match Command::new("nvcc").arg("--list-gpu-arch").output() {
        Ok(o) => o,
        Err(_) => return false,
    };
    if !out.status.success() {
        return false;
    }
    let s = String::from_utf8_lossy(&out.stdout);
    s.lines().any(|l| l.trim() == arch)
}

fn main() {
    let want_cuda = env::var("CARGO_FEATURE_CUDA").is_ok();

    if !want_cuda {
        return;
    }

    if !has_nvcc() {
        println!("cargo:warning=nvcc not found — CUDA disabled");
        return;
    }

    println!("cargo:warning=CUDA enabled");

    println!("cargo:rerun-if-changed=cuda/ntt_kernel.cu");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let cuda_src = "cuda/ntt_kernel.cu";
    let cuda_obj = out_dir.join("ntt_kernel.o");
    let cuda_lib = out_dir.join("libntt_cuda.a");

    let mut compute_caps: Vec<&str> = vec![
        "compute_75", // RTX 20xx
        "compute_86", // RTX 30xx
        "compute_89", // RTX 40xx (Ada)
    ];

    // RTX 50xx (Blackwell, sm_120) — only emit native code when the toolkit
    // supports it. Otherwise the 5090 falls back to JIT-compiling the
    // forward-compat PTX below, which still works but pays a one-time cost.
    let blackwell = nvcc_supports_arch("compute_120");
    if blackwell {
        compute_caps.push("compute_120");
    }

    let mut gencode = Vec::new();
    for cc in &compute_caps {
        let sm = cc.replace("compute_", "sm_");
        gencode.push("-gencode".into());
        gencode.push(format!("arch={},code={}", cc, sm));
    }

    // Forward-compat PTX of the highest available arch so even newer GPUs
    // can JIT it.
    let highest_ptx = if blackwell { "compute_120" } else { "compute_89" };
    gencode.push("-gencode".into());
    gencode.push(format!("arch={0},code={0}", highest_ptx));

    let status = Command::new("nvcc")
        .arg("-c")
        .arg(cuda_src)
        .arg("-o")
        .arg(&cuda_obj)
        .arg("-Xcompiler")
        .arg("-fPIC")
        .arg("-O3")
        .args(&gencode)
        .status()
        .expect("Failed to spawn nvcc");

    if !status.success() {
        panic!("nvcc failed");
    }

    let status = Command::new("ar")
        .arg("rcs")
        .arg(&cuda_lib)
        .arg(&cuda_obj)
        .status()
        .expect("Failed to spawn ar");

    if !status.success() {
        panic!("ar failed");
    }

    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=ntt_cuda");

    let cuda_lib_paths = [
        "/usr/local/cuda/lib64",
        "/usr/lib/x86_64-linux-gnu",
        "/usr/lib64",
    ];

    for path in cuda_lib_paths {
        println!("cargo:rustc-link-search=native={}", path);
    }

    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-cfg=has_cuda");
}
