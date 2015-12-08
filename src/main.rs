extern crate rustc_serialize;
extern crate docopt;
extern crate num;
extern crate flate2;
extern crate byteorder;
extern crate bincode;

mod tallyrun;
mod tallyread;
mod seqtable;

use docopt::Docopt;
use std::io::prelude::*;
use std::fs::File;
use std::io::BufReader;
use std::process::exit;
use std::collections::VecDeque;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use flate2::read::GzDecoder;

enum State {
    HeaderStart,
    HeaderChrom,
    Header,
    Body,
}

#[derive(Debug,Clone)]
struct FreqPair {
    plus: u64,
    minus: u64,
}

#[derive(Debug)]
struct UnMap {
    plus: BTreeSet<u64>,
    minus: BTreeSet<u64>,
}

struct EnzContext<'a> {
    cut_dna_value: u64,
    cut_size: usize,
    
    // auxiliary data
    radix_power: u64,
    buf: VecDeque<u8>,
    
    // filters
    plus_offset: u64,
    minus_offset: u64,
    pluspos: u64,
    minuspos: u64,
    sequnmap: Option<&'a UnMap>,
}

impl<'a> EnzContext<'a> {
    // private API
    fn update_cut_dna_value(&mut self, value: u8) -> bool {
        if self.buf.len() == self.cut_size {
            if let Some(head) = self.buf.pop_front() {
                self.cut_dna_value -= self.radix_power * (head as u64);
            }
        }
        self.buf.push_back(value);
        self.cut_dna_value *= 4;
        self.cut_dna_value += value as u64;
    
        return self.buf.len() == self.cut_size;
    }
    
    fn update_freqs(&mut self, freqs: &mut Vec<FreqPair>) {
        if let Some(umap) = self.sequnmap {
            if ! umap.plus.contains(&self.pluspos) {
                freqs[self.cut_dna_value as usize].plus += 1;
            }
            if ! umap.minus.contains(&self.minuspos) {
                freqs[self.cut_dna_value as usize].minus += 1;
            }
        } else {
            freqs[self.cut_dna_value as usize].plus += 1;
            freqs[self.cut_dna_value as usize].minus += 1;
        }
    }
    
    // public API
    pub fn new(cut_size: usize, plus_offset: u64, minus_offset: u64) -> EnzContext<'a> {
        EnzContext {
            cut_dna_value: 0u64,
            cut_size: cut_size,
            radix_power: num::pow(4, cut_size - 1),
            buf: VecDeque::new(),
            plus_offset: plus_offset,
            minus_offset: minus_offset,
            pluspos: plus_offset,
            minuspos: cut_size as u64 - minus_offset,
            sequnmap: None,
        }
    }
    
    pub fn sequence_change(&mut self, sequnmap: Option<&'a UnMap>) {
        self.pluspos = self.plus_offset;
        self.minuspos = self.cut_size as u64 - self.minus_offset;
        self.sequnmap = sequnmap;
        self.buf.clear();
        self.cut_dna_value = 0;
    }
    
    pub fn add_base(&mut self, base: u8, freqs: &mut Vec<FreqPair>) {
        if self.update_cut_dna_value(base) {
            self.update_freqs(freqs);
        }
        self.pluspos += 1;
        self.minuspos += 1;
    }
    
    pub fn add_n(&mut self) {
        self.pluspos += 1;
        self.minuspos += 1; 
        self.cut_dna_value = 0;
        self.buf.clear();
    }
}

fn get_unmap(unmap: &Option<BTreeMap<usize,UnMap>>, idx: usize) -> Option<&UnMap> {
    if unmap.is_some() {
        unmap.as_ref().unwrap().get(&idx)
    } else {
        None
    }
}

fn process_fasta<R: Read>(reader: R, cut_size: usize, unmap: Option<BTreeMap<usize,UnMap>>, plus_offset: u64, minus_offset: u64) -> Vec<FreqPair> {
    let mut iter = reader.bytes();//.take(70);
    let mut state = State::HeaderStart;
    let mut enzctxt = EnzContext::new(cut_size, plus_offset, minus_offset);
    let mut freqs: Vec<FreqPair> = std::iter::repeat(FreqPair{ plus: 0u64, minus: 0u64 }).take(num::pow(4, cut_size)).collect();
    let mut chrom: Vec<u8> = Vec::new();
    let mut seqidx: isize = -1;
    let mut seqpos = 0u64;
    
    while let Some(Ok(byte)) = iter.next() {
        match state {
            State::HeaderStart => if byte == b'>' {
                    chrom.clear();
                    state = State::HeaderChrom;
                    seqidx += 1;
                    seqpos = 0;
                    enzctxt.sequence_change(get_unmap(&unmap, seqidx as usize));
                } else {
                    println!("Invalid FASTA file. {}", byte as char);
                    exit(1);
                },
            State::HeaderChrom =>
                if byte == b' ' {
                    state = State::Header;
                } else if byte == b'\n' {
                    println!("{:?}", String::from_utf8_lossy(&chrom)); 
                    state = State::Body;
                } else {
                    chrom.push(byte);
                },
            State::Header => if byte == b'\n' { println!("[{}] {:?}", seqidx, String::from_utf8_lossy(&chrom)); state = State::Body },
            State::Body => {
                match byte {
                    b'>' => {
                        println!(" - {} bases", seqpos + 1);
                        chrom.clear();
                        state = State::HeaderChrom;
                        seqidx += 1;
                        seqpos = 0;
                        enzctxt.sequence_change(get_unmap(&unmap, seqidx as usize));
                    },
                    b'a' | b'A' => { enzctxt.add_base(0, &mut freqs); seqpos += 1; },
                    b'c' | b'C' => { enzctxt.add_base(1, &mut freqs); seqpos += 1; },
                    b'g' | b'G' => { enzctxt.add_base(2, &mut freqs); seqpos += 1; },
                    b't' | b'T' => { enzctxt.add_base(3, &mut freqs); seqpos += 1; },
                    b'n' | b'N' => {
                        seqpos += 1;
                        enzctxt.add_n();
                    },
                    _ => {},
               }
            },
        };
    }
    return freqs;
}

fn process_file(path: String, cut_size: usize, unmap: Option<BTreeMap<usize,UnMap>>, plus_offset: u64, minus_offset: u64) -> Vec<FreqPair> {
    println!("hello");
    let f = File::open(path.clone()).ok().expect("Can't open file.");
    
    match GzDecoder::new(f) {
        Ok(reader) => process_fasta(BufReader::new(reader), cut_size, unmap, plus_offset, minus_offset),
        Err(_) => {
            // re-open file
            let f = File::open(path).ok().expect("Can't open file.");
            let reader = BufReader::new(f);
            process_fasta(reader, cut_size, unmap, plus_offset, minus_offset)
        },
    }
}

fn parse_line(bytes: Vec<u8>) -> Result<(usize,  i64), std::io::Error> {
    use std::io::{Error, ErrorKind};
    let mut seq_idx = 0usize;
    let mut pos = 0i64;
    let mut idx = 0;
    
    while bytes[idx] >= b'0' && bytes[idx] <= b'9' {
        seq_idx = seq_idx * 10 + (bytes[idx] - b'0') as usize;
        idx += 1;
    }
    
    if bytes[idx] != b'\t' {
        return Err(Error::new(ErrorKind::Other, "unexpected byte in tallymer line"));
    }
    
    // parse sign
    idx += 1;
    let is_minus = match bytes[idx] {
        b'-' => true,
        b'+' => false,
        _ => return Err(Error::new(ErrorKind::Other, "unexpected byte in tallymer line"))
    };
    idx += 1;
    
    // parse position
    while bytes[idx] >= b'0' && bytes[idx] <= b'9' {
        pos = pos * 10 + (bytes[idx] - b'0') as i64;
        idx += 1;
    }
    
    if is_minus {
        Ok((seq_idx, -pos))
    } else {
        Ok((seq_idx, pos))
    }
}

/// Create one BTreeSet per sequence, allowing efficient query of sparse unmappable set.
fn read_mappability<R: BufRead>(reader: R) -> std::io::Result<BTreeMap<usize,UnMap>> {
    let mut result = BTreeMap::new();
    let mut plus = BTreeSet::new();
    let mut minus = BTreeSet::new();
    let mut prev_seq = Option::None;
    
    // iterate over file
    // expect lines of the form: <number> [+-]<number>
    for line in reader.lines() {
        let line = try!(line);
        
        if let Ok((seq, pos)) = parse_line(line.into_bytes()) {      
            match prev_seq {
                None => {
                    prev_seq = Some(seq);
                },
                Some(val) if val != seq => {
                    result.insert(val, UnMap { plus: plus, minus: minus });
                    //result.push(UnMap { plus: plus, minus: minus });
                    plus = BTreeSet::new();
                    minus = BTreeSet::new();
                    prev_seq = Some(seq);
                },
                _ => {},
            };
            
            // todo: revise to get proper offsets ...
            if pos >= 0 {
                plus.insert(pos as u64);
            } else {
                minus.insert(-pos as u64);
            }
        }
    }
    
    // last entry
    match prev_seq {
        Some(val) => {result.insert(val, UnMap { plus: plus, minus: minus });}, // result.push(UnMap { plus: plus, minus: minus }),
        _ => {},
    };

    return Ok(result);
}

fn load_mappability(path: String) -> std::io::Result<BTreeMap<usize,UnMap>> {
    println!("hello");
    let f = File::open(path.clone()).ok().expect("Can't open file.");
    
    match GzDecoder::new(f) {
        Ok(reader) => read_mappability(BufReader::new(reader)),
        Err(_) => {
            // re-open file
            let f = File::open(path).ok().expect("Can't open file.");
            let reader = BufReader::new(f);
            read_mappability(reader)
        },
    }
}

/* Main usage/arguments */

const USAGE: &'static str = "
Cut-site frequencies

Usage:
  enzcut tallymer <fasta-file> <read-size> [--parts=<n>]
  enzcut <fasta-file> [options]
  enzcut (-h | --help)
  enzcut --version

Options:
  -h --help           Show this screen.
  --version           Show version.
  --cut-size=<n>      Cut-site size [default: 4].
  --tallymer=<file>   Unmappable positions file produced by tallymer (seq, pos).
  --plus-offset=<p>   Cut-site offset on plus strand [default: 2]. Eg, p=2 AA[A]A.
  --minus-offset=<m>  Cut-site offset on minus strand [default: 2]. Eg, m=2 A[A]AA.
  --read-size=<r>     Read length [default: 36].
  --parts=<n>         Split suffix tree generation into n parts [default: 4].
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_fasta_file: String,
    arg_read_size: u16,
    flag_cut_size: usize,
    flag_tallymer: Option<String>,
    flag_plus_offset: u64,
    flag_minus_offset: u64,
    flag_version: bool,
    flag_read_size: u16,
    flag_parts: u8,
    cmd_tallymer: bool,
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
    
    if args.cmd_tallymer {
        let tally_path = tallyrun::tallymer_createfile(&args.arg_fasta_file, args.arg_read_size, args.flag_parts);
        println!("produced/found {:}", tally_path.to_str().unwrap());
        
        return;
    }
    
    let tally_path = tallyrun::tallymer_createfile(&args.arg_fasta_file, args.flag_read_size, args.flag_parts);
    println!("produced/found {:}", tally_path.to_str().unwrap());
    
    
    println!("do more stuff");
    /*
    
    let unmap = if let Some(path) = args.flag_tallymer {
        if let Ok(value) = load_mappability(path) {
            Some(value)
        } else {
            None
        }
    } else {
        None
    };
    
    // process file
    let freqs = process_file(args.arg_fasta_file, args.flag_cut_size, unmap, args.flag_plus_offset, args.flag_minus_offset);

    // output
    println!("Yay! freqs: {:?}", freqs);
    let sum = freqs.iter().fold(0, |sum, x| sum + x.plus);
    println!("Sum: {}", sum);*/
}
