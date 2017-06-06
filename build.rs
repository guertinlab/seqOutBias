// build.rs

extern crate cheddar;
extern crate rustc_version;

use rustc_version::{version, Version};
use std::io::{self, Write};
use std::process::exit;

fn main() {
    cheddar::Cheddar::new().expect("could not read manifest")
        .module("c_api").expect("malformed module path")
        .run_build("target/include/seqoutbiaslib.h");
    
    // Check for a minimum version
    if version().unwrap() < Version::parse("1.11.0").unwrap() {
        writeln!(&mut io::stderr(), "This crate requires rustc >= 1.11.0.").unwrap();
        exit(1);
    }
}
