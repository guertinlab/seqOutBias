//!
//!	Core library entry point.
//!
extern crate rustc_serialize;
extern crate byteorder;
extern crate bincode;
extern crate flate2;
extern crate rust_htslib as htslib;

use std::fs;
use std::io::ErrorKind;
use std::process::exit;

// generic functions
pub fn file_exists<P: AsRef<Path>>(filename: P) -> bool {
	let ref_filename = filename.as_ref();
	match fs::metadata(ref_filename) {
		Ok(meta) => meta.is_file(),
		Err(err) => match err.kind() {
			ErrorKind::NotFound => false,
			_ => { println!("Failed to open file {}: {:}", ref_filename.to_string_lossy(), err); exit(-1) },
		}, 
	}
}

// private modules
mod randfile;

// public modules
pub mod tallyrun;
pub mod tallyread;
pub mod seqtable;
pub mod fasta;
pub mod filter;
pub mod counts;
pub mod bigwig;
pub mod scale;
pub mod outputfile;

// C API
pub mod c_api;
pub use c_api::*;
use std::path::Path;
