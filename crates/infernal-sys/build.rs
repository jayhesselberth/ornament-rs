//! Build script for compiling Infernal, HMMER, and Easel C libraries

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    // Path to ext/ directory (two levels up from crates/infernal-sys)
    let ext_dir = manifest_dir.parent().unwrap().parent().unwrap().join("ext");
    let _ = manifest_dir; // unused now
    let infernal_dir = ext_dir.join("infernal");

    // HMMER and Easel are inside the infernal directory
    let easel_dir = infernal_dir.join("easel");
    let hmmer_dir = infernal_dir.join("hmmer");

    // Check that dependencies are set up
    check_dependencies(&infernal_dir, &hmmer_dir, &easel_dir);

    // Build output directory - build everything in infernal's build dir
    let build_dir = out_dir.join("infernal-build");
    std::fs::create_dir_all(&build_dir).unwrap();

    // Detect architecture for SIMD flags
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();

    // Rebuild if source changes
    println!("cargo:rerun-if-changed={}", infernal_dir.display());

    // Configure and build Infernal (which builds HMMER and Easel too)
    build_infernal(&infernal_dir, &build_dir, &target_arch);

    // Link the libraries
    println!("cargo:rustc-link-search=native={}", build_dir.join("easel").display());
    println!("cargo:rustc-link-search=native={}", build_dir.join("hmmer/src").display());
    println!("cargo:rustc-link-search=native={}", build_dir.join("src").display());

    println!("cargo:rustc-link-lib=static=infernal");
    println!("cargo:rustc-link-lib=static=hmmer");
    println!("cargo:rustc-link-lib=static=easel");

    // System libraries
    println!("cargo:rustc-link-lib=m");

    // Generate bindings
    generate_bindings(&infernal_dir, &hmmer_dir, &easel_dir, &out_dir);
}

fn check_dependencies(infernal_dir: &Path, hmmer_dir: &Path, easel_dir: &Path) {
    let mut missing = Vec::new();

    if !infernal_dir.exists() {
        missing.push("ext/infernal");
    }
    if !hmmer_dir.exists() {
        missing.push("ext/infernal/hmmer");
    }
    if !easel_dir.exists() {
        missing.push("ext/infernal/easel");
    }

    if !missing.is_empty() {
        eprintln!("\n");
        eprintln!("=======================================================");
        eprintln!("ERROR: Missing dependencies!");
        eprintln!("=======================================================");
        eprintln!("");
        eprintln!("The following directories are missing:");
        for m in &missing {
            eprintln!("  - {}", m);
        }
        eprintln!("");
        eprintln!("Please run the setup script to clone dependencies:");
        eprintln!("");
        eprintln!("  ./scripts/setup-deps.sh");
        eprintln!("");
        eprintln!("Or manually clone:");
        eprintln!("  git clone -b develop https://github.com/EddyRivasLab/infernal.git ext/infernal");
        eprintln!("  cd ext/infernal");
        eprintln!("  git clone -b develop https://github.com/EddyRivasLab/hmmer.git");
        eprintln!("  git clone -b develop https://github.com/EddyRivasLab/easel.git");
        eprintln!("  autoconf");
        eprintln!("");
        eprintln!("=======================================================\n");
        panic!("Missing dependencies. See instructions above.");
    }

    // Check for configure script
    if !infernal_dir.join("configure").exists() {
        eprintln!("\n");
        eprintln!("=======================================================");
        eprintln!("ERROR: Infernal not configured!");
        eprintln!("=======================================================");
        eprintln!("");
        eprintln!("Please run autoconf in the infernal directory:");
        eprintln!("");
        eprintln!("  cd ext/infernal && autoconf");
        eprintln!("");
        eprintln!("Or run the setup script:");
        eprintln!("");
        eprintln!("  ./scripts/setup-deps.sh");
        eprintln!("");
        eprintln!("=======================================================\n");
        panic!("Infernal not configured. See instructions above.");
    }
}

fn build_infernal(src_dir: &Path, build_dir: &Path, target_arch: &str) {
    // Check if already built
    let lib_path = build_dir.join("src/libinfernal.a");
    if lib_path.exists() {
        println!("cargo:warning=Infernal already built, skipping...");
        return;
    }

    // Configure
    let mut configure_cmd = Command::new(src_dir.join("configure"));
    configure_cmd
        .current_dir(build_dir)
        .arg(format!("--prefix={}", build_dir.display()));

    // Add SIMD flags based on architecture
    match target_arch {
        "x86_64" => {
            configure_cmd.arg("--enable-sse");
        }
        "aarch64" => {
            // NEON is auto-detected on ARM64
        }
        _ => {}
    }

    run_command(&mut configure_cmd, "configure infernal");

    // Make (builds Infernal, HMMER, and Easel)
    run_command(
        Command::new("make")
            .current_dir(build_dir)
            .arg("-j")
            .arg(num_cpus().to_string()),
        "make infernal"
    );
}

fn generate_bindings(infernal_dir: &Path, hmmer_dir: &Path, easel_dir: &Path, out_dir: &Path) {
    // Build directory contains generated config headers (esl_config.h, etc.)
    let build_dir = out_dir.join("infernal-build");

    let bindings = bindgen::Builder::default()
        // Infernal main header
        .header(infernal_dir.join("src/infernal.h").to_string_lossy())
        // Include paths - source directories
        .clang_arg(format!("-I{}", infernal_dir.join("src").display()))
        .clang_arg(format!("-I{}", hmmer_dir.join("src").display()))
        .clang_arg(format!("-I{}", easel_dir.display()))
        // Include paths - build directories (for generated config headers)
        .clang_arg(format!("-I{}", build_dir.join("src").display()))
        .clang_arg(format!("-I{}", build_dir.join("hmmer/src").display()))
        .clang_arg(format!("-I{}", build_dir.join("easel").display()))
        // Allowlist the types/functions we need
        // Infernal types
        .allowlist_type("CM_t")
        .allowlist_type("CM_FILE")
        .allowlist_type("CM_TOPHITS")
        .allowlist_type("CM_HIT")
        .allowlist_type("CM_PIPELINE")
        .allowlist_type("Parsetree_t")
        // HMMER types (for HMM filter in pipeline)
        .allowlist_type("P7_OPROFILE")
        .allowlist_type("P7_PROFILE")
        .allowlist_type("P7_BG")
        .allowlist_type("P7_HMM")
        .allowlist_type("P7_SCOREDATA")
        // Easel types
        .allowlist_type("ESL_ALPHABET")
        .allowlist_type("ESL_SQ")
        .allowlist_type("ESL_SQFILE")
        .allowlist_type("ESL_DSQ")
        .allowlist_type("ESL_MSA")
        .allowlist_type("ESL_GETOPTS")
        // Infernal functions
        .allowlist_function("cm_file_Open")
        .allowlist_function("cm_file_Read")
        .allowlist_function("cm_file_Close")
        .allowlist_function("cm_Configure")
        .allowlist_function("cm_pipeline_Create")
        .allowlist_function("cm_pipeline_Destroy")
        .allowlist_function("cm_Pipeline")
        .allowlist_function("cm_tophits_.*")
        .allowlist_function("CreateCMConsensus")
        .allowlist_function("FreeCM")
        // HMMER functions (for HMM filter setup)
        .allowlist_function("p7_profile_Create")
        .allowlist_function("p7_profile_Destroy")
        .allowlist_function("p7_oprofile_Create")
        .allowlist_function("p7_oprofile_Convert")
        .allowlist_function("p7_oprofile_Destroy")
        .allowlist_function("p7_bg_Create")
        .allowlist_function("p7_bg_Destroy")
        .allowlist_function("p7_ProfileConfig")
        .allowlist_function("p7_hmm_ScoreDataCreate")
        .allowlist_function("p7_hmm_ScoreDataDestroy")
        // Easel functions
        .allowlist_function("esl_alphabet_Create")
        .allowlist_function("esl_alphabet_Destroy")
        .allowlist_function("esl_sq_.*")
        .allowlist_function("esl_sqfile_.*")
        .allowlist_function("esl_sqio_.*")
        .allowlist_function("esl_vec_FCopy")
        // Generate
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

fn run_command(cmd: &mut Command, description: &str) {
    println!("cargo:warning=Running: {}", description);

    let status = cmd
        .status()
        .unwrap_or_else(|e| panic!("Failed to execute {}: {}", description, e));

    if !status.success() {
        panic!("{} failed with status: {}", description, status);
    }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1)
}
