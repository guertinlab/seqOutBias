// build.rs
#[cfg(feature = "headers")]
extern crate cbindgen;
extern crate rustc_version;

use rustc_version::version_matches;   
use std::io::{self, Write};
use std::process::exit;
#[cfg(feature = "headers")]
use std::env;
#[cfg(feature = "headers")]
use std::path::PathBuf;

fn main() {

    #[cfg(feature = "headers")]
    {
        let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
        let package_name = "seqoutbiaslib";

        let output_file = target_dir()
        .join( "include" )
        .join(format!("{}.h", package_name))
        .display()
        .to_string();

        let config = cbindgen::Config {
          language: cbindgen::Language::C,
          include_guard: Some(format!("cbindgen_generated_{}_h", package_name)),
          cpp_compat: true,
          style: cbindgen::Style::Type,
          ..Default::default()
        };

        cbindgen::generate_with_config(&crate_dir, config)
        .unwrap()
        .write_to_file(&output_file);
    }

    // Check for a minimum version
    if !version_matches(">= 1.32.0") {
        writeln!(&mut io::stderr(), "This crate requires rustc >= 1.32.0.").unwrap();
        exit(1);
    }
}

#[cfg(feature = "headers")]
fn target_dir() -> PathBuf {
    if let Ok(target) = env::var("CARGO_TARGET_DIR") {
        PathBuf::from(target)
    } else {
        PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap()).join("target")
    }
}
