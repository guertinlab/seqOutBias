// build.rs

extern crate cheddar;

fn main() {
    cheddar::Cheddar::new().expect("could not read manifest")
        .module("c_api").expect("malformed module path")
        .run_build("target/include/seqoutbiaslib.h");
}
