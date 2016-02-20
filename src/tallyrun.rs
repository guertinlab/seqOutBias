//!
//!	This module is responsible for running tallymer to create a mappability file
//!	if no such file is present.
//!
extern crate flate2;

use std::process::Command;
use std::process::exit;
use std::process::Stdio;
use std::process::Child;
use std::path::Path;
use std::ffi::OsString;
use std::ffi::OsStr;
use std::fs;
use std::fs::File;
use std::io::Result;
use std::io::ErrorKind;
use std::io::copy;
use std::io::BufWriter;
use flate2::write::GzEncoder;
use flate2::Compression;


/// Get basename for file after excluding a .gz suffix if it exists
fn basename_nogz(filepath: &str) -> OsString {
	let filepath = if filepath.ends_with(".gz") {
		&filepath[0..(filepath.len() - 3)]
	} else {
		filepath
	};
	
	let path = Path::new(filepath);
	path.file_stem().unwrap_or(path.file_name().expect("Failed to extract filename from path.")).to_os_string()
}


/// Create filename with format: <fasta basename>.tal_<readlen>.gtTxt
fn create_output_filename(fasta: &str, readlen: u16) -> OsString {
	let mut stem = basename_nogz(fasta);
	stem.push(format!(".tal_{}.gtTxt", readlen));
	
	return stem;
}

/// Create filename with format: <fasta basename>.sft
fn create_tallymerindex_filename(fasta: &str, readlen: u16) -> OsString {
	let mut stem = basename_nogz(fasta);
	stem.push(format!(".tal_{}", readlen));
	
	return stem;
}


/// Create filename with format: <fasta basename>.sft
fn create_suffixtree_filename(fasta: &str) -> OsString {
	let mut stem = basename_nogz(fasta);
	stem.push(".sft");
	
	return stem;
}

/// Check if file exists by reading it's metadata
fn file_exists(filename: &OsStr) -> bool {
	match fs::metadata(filename) {
		Ok(meta) => meta.is_file(),
		Err(err) => match err.kind() {
			ErrorKind::NotFound => false,
			_ => { println!("Failed to open file {}: {:}", filename.to_str().unwrap(), err); exit(-1) },
		}, 
	}
}

fn wait_or_exit(spawned: Result<Child>, name: &str) {
	match spawned {
		Ok(mut child) => {
			let ecode = child.wait().unwrap_or_else(|e| {
				println!("failed to wait on '{}': {}", name, e);
				exit(-1)
			});
			if !ecode.success() {
				println!("'{}' failed!\n", name);
				exit(-1);
			}
		},
		Err(e) => {
			println!("failed to execute '{}': {}", name, e);
			exit(-1)
		}
	}
}


fn pipe_to_gzipped_file_or_exit(proc1: Result<Child>, output: &OsStr) {
	let mut proc1 = proc1.ok().expect("spawn gt tallymer");
	
	let file = File::create(output).ok().expect("create file");
	let writer = BufWriter::new(file);
	let mut encoder = GzEncoder::new(writer, Compression::Best);
	copy(proc1.stdout.as_mut().unwrap(), &mut encoder).ok().expect("write to file");
	
	let ecode = proc1.wait()
    	.unwrap_or_else(|e| {
			println!("Failed to wait on 'gt tallymer': {}", e);
			exit(-1)
		});

	if !ecode.success() {
			println!("gt tallymer!\n");
			exit(-1);
	}
}

/// Create tallymer based mapabillity file 
fn tallymer_create(fasta: &str, readlen: u16, parts: u8, outfile: &OsStr) {
	let sft_filename = create_suffixtree_filename(fasta);
	let tidx_filename = create_tallymerindex_filename(fasta, readlen);
	
	// create suffix-tree index if it does not exist
	let mut tmp = sft_filename.clone();
	tmp.push(".suf"); // actually there are a bunch of files that must exist ...
	if !file_exists(&tmp) {
		let gt_proc = Command::new("gt")
			.arg("suffixerator")
			.arg("-dna")
			.arg("-pl")
			.arg("-tis")
			.arg("-suf")
			.arg("-lcp")
			.arg("-v")
			.arg("-parts")
			.arg(format!("{}", parts)) // reduce in-memory requirements to 1/nth
			.arg("-db")
			.arg(fasta)
			.arg("-indexname")
			.arg(&sft_filename)
			.spawn();
		
		wait_or_exit(gt_proc, "gt suffixerator");
	}
	
	// create genome index
	let gt_proc = Command::new("gt")
		.arg("tallymer")
		.arg("mkindex")
		.arg("-mersize")
		.arg(format!("{}", readlen))
		.arg("-minocc")
		.arg("2")
		.arg("-indexname")
		.arg(&tidx_filename)
		.arg("-counts")
		.arg("-pl")
		.arg("-esa")
		.arg(&sft_filename)
		.spawn();
	
	wait_or_exit(gt_proc, "gt tallymer mkindex");
	
	// generate final file
	let gt_proc = Command::new("gt")
		.arg("tallymer")
		.arg("search")
		.arg("-output")
		.arg("qseqnum")
		.arg("qpos")
		.arg("-strand")
		.arg("fp")
		.arg("-tyr")
		.arg(&tidx_filename)
		.arg("-q")
		.arg(fasta)
		.stdout(Stdio::piped())
		.spawn();
	
	pipe_to_gzipped_file_or_exit(gt_proc, outfile);
}

/// Get tallymer mappability file's path
///
/// This will invoke the 'tallymer' program from genome tools to create the
/// file if it's not found.
///
/// Filename: <fasta basename>.tal_<readlen>.gtTxt
/// Files can optionally have a .gz extension.
///
/// The generated file will be compressed with gzip if the input fasta has an
/// .gz extension.
pub fn tallymer_createfile(fasta: &str, readlen: u16, parts: u8) -> OsString {
	
	// try to locate file
	let mut output = create_output_filename(fasta, readlen);
	
	if file_exists(&output) {
		return output;
	}
	
	// try with .gz extension
	output.push(".gz");
	if file_exists(&output) {
		return output;
	}
	
	// not found, let's create one
	println!("################### Creating mappability file using tallymer ###################");
	tallymer_create(fasta, readlen, parts, &output);
	println!("################################################################################");
	
	return output;
}
