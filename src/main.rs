extern crate rustc_serialize;
extern crate docopt;
extern crate num;

use docopt::Docopt;
use std::io::prelude::*;
use std::fs::File;
use std::io::BufReader;
use std::collections::VecDeque;

enum State {
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

fn process_file<R: Read>(reader: R, cut_size: usize) -> Vec<u64> {
    let mut iter = reader.bytes();//.take(70);
    let mut state = State::Header;
    let mut buf = VecDeque::new();
    let mut tetra = 0u64;
    let radix_power = num::pow(4, cut_size - 1);
    let mut freqs: Vec<u64> = std::iter::repeat(0u64).take(num::pow(4, cut_size)).collect();
    
    while let Some(Ok(byte)) = iter.next() {
        match state {
            State::Header => if byte == b'\n' { state = State::Body },
            State::Body => {
                match byte {
                    b'>' => {
                        state = State::Header;
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
    
    // load data
    let f = File::open(args.arg_fasta_file).ok().expect("Can't open file.");
    let reader = BufReader::new(f);
    
    // process file
    let freqs = process_file(reader, args.flag_cut_size);

    // output
    println!("Yay! freqs: {:?}", freqs);
    let sum = freqs.iter().fold(0, |sum, x| sum + x);
    println!("Sum: {}", sum);
}
