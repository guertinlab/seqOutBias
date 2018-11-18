//!
//!	This module has utility functions to generate random temporary filenames
//!	It replaces tempfile usage which does not seem to be working out as intended.
//!
extern crate rand;
use std::path::{PathBuf};

use self::rand::{thread_rng, Rng};

fn random_string(size: usize) -> String {
    let rand_string: String = thread_rng()
        .gen_ascii_chars()
        .take(size)
        .collect();
    rand_string
}

pub fn random_file(prefix: &str, suffix: &str, rand_bytes: usize) -> PathBuf {
    loop {
        let rand_string = random_string(rand_bytes);
        let mut filename = prefix.to_owned();
        filename.push_str(&rand_string);
        filename.push_str(suffix);

        let mut path = PathBuf::new();
        path.push(filename);
        if !path.exists() {
            return path;
        }
    }
}
