//!
//! Code to read FASTA file
//!
use std::collections::VecDeque;
use std::io::prelude::*;
use seqtable::{SeqTableParams,SeqTableWriter,SequenceWriter,SeqBuffer};
use tallyread::UnMap;
use std::fs::File;
use std::io::{BufReader,Bytes};
use std::ffi::OsStr;
use flate2::read::GzDecoder;
use std::process::exit;

/// Holds information on the current enzyme cut region
struct EnzContext {
    cut_dna_value: u64,
    cut_size: usize,
    
    // auxiliary data
    radix_power: u64,
    buf: VecDeque<u8>,
}

impl EnzContext {
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
        
    // public API
    pub fn new(cut_size: u8) -> EnzContext {
        EnzContext {
            cut_dna_value: 0u64,
            cut_size: cut_size as usize,
            radix_power: 4u64.pow((cut_size - 1) as u32),
            buf: VecDeque::new(),
        }
    }
    
    pub fn sequence_change(&mut self) {
        self.buf.clear();
        self.cut_dna_value = 0;
    }
    
    pub fn add_base(&mut self, base: u8) -> Option<u64> {
        if self.update_cut_dna_value(base) {
            Some(self.cut_dna_value)
        } else {
            None
        }
    }
    
    pub fn add_n(&mut self) {
        self.cut_dna_value = 0;
        self.buf.clear();
    }
}

// FASTA State machine
#[derive(PartialEq)]
enum State {
    HeaderStart,
    HeaderChrom,
    Header,
    End,
}

macro_rules! store_base {
    ($e:expr, $b:expr, $n:expr) => {
        match $e.add_base($n) {
            Some(nmer_index) => $b.push((nmer_index + 1) as u32),
            None => $b.push(0),
    }};
}

fn process_sequence<R1: Read, R2: BufRead>(seqwrt: SequenceWriter<File>, iter: &mut Bytes<R1>, enzctxt: &mut EnzContext, params: &SeqTableParams, unmap: &UnMap<R2>) -> State {
    let mut buf = SeqBuffer::new(seqwrt, *params, unmap);
    let mut seqpos = 0u32;
    
    while let Some(Ok(byte)) = iter.next() {
        match byte {
            b'>' => {
                println!("# - {} bases", seqpos + 1);                
                return State::HeaderChrom;
            },
            b'a' | b'A' => { store_base!(enzctxt, buf, 0); seqpos += 1; }, 
            b'c' | b'C' => { store_base!(enzctxt, buf, 1); seqpos += 1; },
            b'g' | b'G' => { store_base!(enzctxt, buf, 2); seqpos += 1; },
            b't' | b'T' => { store_base!(enzctxt, buf, 3); seqpos += 1; },
            b'n' | b'N' => {
                enzctxt.add_n();
                buf.push(0);
                seqpos += 1;
            },
            _ => {},
        }
    }
    println!("# - {} bases", seqpos + 1);
    State::End
}

/// Read FASTA file and produce SeqTable file
pub fn generate_seqtable<R1: Read, R2: BufRead>(fasta: R1, tallymer: R2, params: SeqTableParams, outfile: &str) {
	let mut unmap = UnMap::open(tallymer).ok().expect("Load mappability information");
    let mut iter = fasta.bytes();
    let mut state = State::HeaderStart;
    let mut enzctxt = EnzContext::new(params.cut_length);
    let mut chrom: Vec<u8> = Vec::new();
    
    let f_out = File::create(outfile).ok().expect("create file");
    let mut output = SeqTableWriter::new(f_out, params, 3200000).ok().expect("create store");
    
    while let Some(Ok(byte)) = iter.next() {
        match state {
            State::HeaderStart => if byte == b'>' {
                    chrom.clear();
                    state = State::HeaderChrom;
                    enzctxt.sequence_change();
                } else {
                    println!("Invalid FASTA file. {}", byte as char);
                    exit(1);
                },
            State::HeaderChrom =>
                if byte == b' ' {
                    state = State::Header;
                } else if byte == b'\n' {
                    let seqwrt = output.create_sequence(String::from_utf8_lossy(&chrom).into_owned());
                    println!("# chrom: {:?}", String::from_utf8_lossy(&chrom)); 
                    state = process_sequence(seqwrt, &mut iter, &mut enzctxt, &params, &mut unmap);
                    
                    // after processing sequence
                    if state == State::HeaderChrom {
                        chrom.clear();
                        enzctxt.sequence_change();
                        unmap.read_next_sequence().ok().expect("failed to read tallymer data");
                    }
                } else {
                    chrom.push(byte);
                },
            State::Header => if byte == b'\n' {
                    let seqwrt = output.create_sequence(String::from_utf8_lossy(&chrom).into_owned());
                    println!("# chrom: {:?}", String::from_utf8_lossy(&chrom)); 
                    state = process_sequence(seqwrt, &mut iter, &mut enzctxt, &params, &mut unmap);
                    
                    // after processing sequence
                    if state == State::HeaderChrom {
                        chrom.clear();
                        enzctxt.sequence_change();
                        unmap.read_next_sequence().ok().expect("failed to read tallymer data");
                    }
                },
            State::End => {},
        };
    }
}

pub fn process_fasta(fasta_path: &str, tallymer_path: &OsStr, params: SeqTableParams, outfile: &str) {
    let f_fasta = File::open(fasta_path).ok().expect("Can't open FASTA file.");
    let f_tallymer = File::open(tallymer_path).ok().expect("Can't open Tallymer file.");
    
    match GzDecoder::new(f_fasta) {
        Ok(reader_fasta) => {
            match GzDecoder::new(f_tallymer) {
                Ok(reader_tallymer) => generate_seqtable(reader_fasta, BufReader::new(reader_tallymer), params, outfile),
                Err(_) => {
                    // re-open file
                    let f_tallymer = File::open(tallymer_path).ok().expect("Can't open Tallymer file.");
                    let reader_tallymer = BufReader::new(f_tallymer);
                    
                    generate_seqtable(reader_fasta, reader_tallymer, params, outfile);
                },
            }
        },
        Err(_) => {
            // re-open file
            let f_fasta = File::open(fasta_path).ok().expect("Can't open FASTA file.");
            let reader_fasta = BufReader::new(f_fasta);
            match GzDecoder::new(f_tallymer) {
                Ok(reader_tallymer) => generate_seqtable(reader_fasta, BufReader::new(reader_tallymer), params, outfile),
                Err(_) => {
                    // re-open file
                    let f_tallymer = File::open(tallymer_path).ok().expect("Can't open Tallymer file.");
                    let reader_tallymer = BufReader::new(f_tallymer);
                    
                    generate_seqtable(reader_fasta, reader_tallymer, params, outfile);
                },
            }
        },
    }
}