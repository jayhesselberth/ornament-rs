//! Raw FFI bindings to Infernal, HMMER, and Easel
//!
//! This crate provides low-level unsafe bindings to the Infernal covariance model
//! library and its dependencies (HMMER and Easel).
//!
//! For safe Rust wrappers, use the `modtrna-core` crate instead.

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(clippy::all)]

// Include the generated bindings
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bindings_exist() {
        // Just verify that the types exist
        // Actual functionality tests will be in modtrna-core
        let _size = std::mem::size_of::<CM_t>();
    }
}
