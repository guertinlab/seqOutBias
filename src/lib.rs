//!
//!	Core library entry point.
//!
extern crate rustc_serialize;
extern crate byteorder;
extern crate bincode;
extern crate flate2;
extern crate rust_htslib as htslib;

pub mod tallyrun;
mod tallyread;
pub mod seqtable;
pub mod fasta;
mod filter;
pub mod counts;
mod bigwig;
pub mod scale;
