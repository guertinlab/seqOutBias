// build.rs
#[cfg(feature = "headers")]
extern crate cheddar;
extern crate rustc_version;

use rustc_version::version_matches;   
use std::io::{self, Write};
use std::process::exit;

fn main() {
    #[cfg(feature = "headers")]
    cheddar::Cheddar::new().expect("could not read manifest")
        .module("c_api").expect("malformed module path")
        .run_build("target/include/seqoutbiaslib.h");
    
    // Check for a minimum version
    if !version_matches(">= 1.11.0") {
        writeln!(&mut io::stderr(), "This crate requires rustc >= 1.11.0.").unwrap();
        exit(1);
    }
}
