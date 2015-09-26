extern crate rustc_serialize;
extern crate docopt;
extern crate num;
#[cfg(feature = "gzip")]
extern crate flate2;

use docopt::Docopt;
use std::io::prelude::*;
use std::fs::File;
use std::io::BufReader;
use std::process::exit;
use std::collections::VecDeque;
#[cfg(feature = "gzip")]
use flate2::read::GzDecoder;

enum State {
    HeaderStart,
    HeaderChrom,
    Header,
    Body,
}

fn update(freqs: &mut Vec<u64>, buf: &mut VecDeque<u8>, tetra: &mut u64, value: u8, cap: usize, radix_power: u64) {
    if buf.len() == cap {
        if let Some(head) = buf.pop_front() {
            *tetra = *tetra - radix_power*(head as u64);
        }
    }
    buf.push_back(value);
    *tetra = *tetra * 4;
    *tetra = *tetra + value as u64;

    if buf.len() == cap {
        freqs[*tetra as usize] += 1;
    }
}

fn process_fasta<R: Read>(reader: R, cut_size: usize) -> Vec<u64> {
    let mut iter = reader.bytes();//.take(70);
    let mut state = State::HeaderStart;
    let mut buf = VecDeque::new();
    let mut tetra = 0u64;
    let radix_power = num::pow(4, cut_size - 1);
    let mut freqs: Vec<u64> = std::iter::repeat(0u64).take(num::pow(4, cut_size)).collect();
    let mut chrom: Vec<u8> = Vec::new();
    
    while let Some(Ok(byte)) = iter.next() {
        match state {
            State::HeaderStart => if byte == b'>' {
                    chrom.clear();
                    state = State::HeaderChrom;
                } else {
                    println!("Invalid FASTA file. {}", byte as char);
                    exit(1);
                },
            State::HeaderChrom =>
                if byte == b' ' {
                    state = State::Header;
                } else {
                    chrom.push(byte);
                },
            State::Header => if byte == b'\n' { println!("{:?}", String::from_utf8_lossy(&chrom)); state = State::Body },
            State::Body => {
                match byte {
                    b'>' => {
                        chrom.clear();
                        state = State::HeaderChrom;
                        tetra = 0;
                        buf.clear();
                    },
                    b'a' | b'A' => update(&mut freqs, &mut buf, &mut tetra, 0, cut_size, radix_power),
                    b'c' | b'C' => update(&mut freqs, &mut buf, &mut tetra, 1, cut_size, radix_power),
                    b'g' | b'G' => update(&mut freqs, &mut buf, &mut tetra, 2, cut_size, radix_power),
                    b't' | b'T' => update(&mut freqs, &mut buf, &mut tetra, 3, cut_size, radix_power),
                    b'n' | b'N' => {
                        tetra = 0;
                        buf.clear();
                    },
                    _ => {},
               }
            },
        };
    }
    return freqs;
}

#[cfg(feature = "gzip")]
fn process_file(path: String, cut_size: usize) -> Vec<u64> {
    println!("hello");
    let f = File::open(path.clone()).ok().expect("Can't open file.");
    
    match GzDecoder::new(f) {
        Ok(reader) => process_fasta(BufReader::new(reader), cut_size),
        Err(_) => {
            // re-open file
            let f = File::open(path).ok().expect("Can't open file.");
            let reader = BufReader::new(f);
            process_fasta(reader, cut_size)
        },
    }
}

#[cfg(not(feature = "gzip"))]
fn process_file(path: String, cut_size: usize) -> Vec<u64> {
    let f = File::open(path).ok().expect("Can't open file.");
    
    let reader = BufReader::new(f);
    process_fasta(reader, cut_size)
}

/* Main usage/arguments */

const USAGE: &'static str = "
Cut-site frequencies

Usage:
  treta <fasta-file> [--cut-size=<n>]
  treta (-h | --help)
  treta --version

Options:
  -h --help       Show this screen.
  --version       Show version.
  --cut-size=<n>  Cut-site size [default: 4].
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_fasta_file: String,
    flag_cut_size: usize,
    flag_version: bool,
}

fn main() {
    // Parse command line arguments
    let args: Args = Docopt::new(USAGE)
                            .and_then(|d| d.decode())
                            .unwrap_or_else(|e| e.exit());

    if args.flag_version {
        println!("Cut-site frequencies, v{}.{}.{}", 
            env!( "CARGO_PKG_VERSION_MAJOR" ),
            env!( "CARGO_PKG_VERSION_MINOR" ),
            env!( "CARGO_PKG_VERSION_PATCH" ) );
        return;
    }
    
    // process file
    let freqs = process_file(args.arg_fasta_file, args.flag_cut_size);

    // output
    println!("Yay! freqs: {:?}", freqs);
    let sum = freqs.iter().fold(0, |sum, x| sum + x);
    println!("Sum: {}", sum);
}
