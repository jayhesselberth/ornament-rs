//! Compile the CUDA kernels with `nvcc` into a static lib and link the CUDA runtime — but only
//! when the `cuda` feature is enabled. Without the feature this is a no-op, so the crate builds
//! on hosts with no CUDA toolchain (login nodes).
//!
//! Overrides via env:
//!   NVCC        path to nvcc            (default: `nvcc` on PATH)
//!   CUDA_ARCH   gencode virtual/real    (default: `sm_80`, the cluster's A30 / Ampere)
//!   CUDA_HOME   CUDA install prefix     (used to find libcudart if not on the default path)

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=kernels/msv.cu");
    println!("cargo:rerun-if-env-changed=NVCC");
    println!("cargo:rerun-if-env-changed=CUDA_ARCH");
    println!("cargo:rerun-if-env-changed=CUDA_HOME");

    if env::var_os("CARGO_FEATURE_CUDA").is_none() {
        return; // feature off: nothing to build, no CUDA required.
    }

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let nvcc = env::var("NVCC").unwrap_or_else(|_| "nvcc".to_string());
    let arch = env::var("CUDA_ARCH").unwrap_or_else(|_| "sm_80".to_string());
    let lib = out_dir.join("libornament_gpu_kernels.a");

    // nvcc --lib builds a static library directly from the .cu (device + host wrappers).
    let status = Command::new(&nvcc)
        .args([
            "--lib",
            "-O3",
            "-std=c++14",
            "-Xcompiler",
            "-fPIC",
            &format!("-arch={arch}"),
            "-o",
        ])
        .arg(&lib)
        .arg("kernels/msv.cu")
        .status()
        .unwrap_or_else(|e| {
            panic!(
                "failed to run nvcc ({nvcc}): {e}\n\
                 The `cuda` feature needs the CUDA toolkit. On this cluster nvcc lives on the GPU \
                 nodes (compgpu01-03); build there, or set NVCC=/path/to/nvcc."
            )
        });
    assert!(status.success(), "nvcc failed to compile kernels/msv.cu");

    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=ornament_gpu_kernels");

    // Link the CUDA runtime + C++ stdlib the wrappers pull in. Cover both the classic `lib64`
    // layout and CUDA 12/13's `targets/<triple>/lib`.
    if let Ok(home) = env::var("CUDA_HOME") {
        println!("cargo:rustc-link-search=native={home}/lib64");
        println!("cargo:rustc-link-search=native={home}/targets/x86_64-linux/lib");
    }
    for dir in [
        "/usr/local/cuda/lib64",
        "/usr/local/cuda/targets/x86_64-linux/lib",
        "/usr/lib64",
        "/usr/lib/x86_64-linux-gnu",
    ] {
        println!("cargo:rustc-link-search=native={dir}");
    }
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}
