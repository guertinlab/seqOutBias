//!
//! Code to read FASTA file
//!
use std::collections::VecDeque;
use std::io::prelude::*;
use seqtable::{SeqBuffer, SeqTableParams, SeqTableWriter, SequenceWriter};
use tallyread::UnMap;
use std::fs::File;
use std::io::{BufReader, Bytes};
use std::ffi::OsStr;
use flate2::read::GzDecoder;
use std::process::exit;
use fasta::context::{KmerIndex, EnzContext, EnzContextMasked, EnzContextMaskedStrandSpecific, EnzContextSimple};

mod context;

// Base encoding
#[repr(u8)]
#[derive(PartialEq)]
pub enum DNABases {
    A,
    C,
    G,
    T,
    N
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
          KmerIndex{ plus: Some(plus), minus: Some(minus) } => $b.push( plus + 1, minus + 1 ),
          KmerIndex{ plus: Some(plus), minus: None } => $b.push( plus + 1, 0 ),
          KmerIndex{ plus: None, minus: Some(minus) } => $b.push( 0, minus + 1 ),
          KmerIndex{ plus: None, minus: None } => $b.push( 0, 0 ),
    }};
}

fn process_sequence<R1: Read, R2: BufRead, T:EnzContext>(seqwrt: SequenceWriter<File>, iter: &mut Bytes<R1>, enzctxt: &mut T, params: &SeqTableParams, unmap: &UnMap<R2>) -> State {
    let mut buf = SeqBuffer::new(seqwrt, params, unmap);
    let mut seqpos = 0u32;

    while let Some(Ok(byte)) = iter.next() {
        match byte {
            b'>' => {
                println!("# - {} bases", seqpos + 1);                
                return State::HeaderChrom;
            },
            b'a' | b'A' => { store_base!(enzctxt, buf, DNABases::A as u8); seqpos += 1; },
            b'c' | b'C' => { store_base!(enzctxt, buf, DNABases::C as u8); seqpos += 1; },
            b'g' | b'G' => { store_base!(enzctxt, buf, DNABases::G as u8); seqpos += 1; },
            b't' | b'T' => { store_base!(enzctxt, buf, DNABases::T as u8); seqpos += 1; },
            b'n' | b'N' => { store_base!(enzctxt, buf, DNABases::N as u8); seqpos += 1; },
            _ => {},
        }
    }
    println!("# - {} bases", seqpos + 1);
    State::End
}

pub fn reverse_complement(mer: u32, kmersize: u8) -> u32 {
    if mer == 0 {
        return 0;
    }
    return 4u32.pow( kmersize as u32 ) + 1 - mer;
}

/// Read FASTA file and produce SeqTable file
fn generate_seqtable_ctxt<R1: Read, R2: BufRead, T: EnzContext>(fasta: R1, tallymer: R2, params: &SeqTableParams, mut enzctxt: T, outfile: &str) {
	let mut unmap = UnMap::open(tallymer).ok().expect("Load mappability information");
    let mut iter = fasta.bytes();
    let mut state = State::HeaderStart;
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
                    state = process_sequence(seqwrt, &mut iter, &mut enzctxt, params, &mut unmap);
                    
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
                    state = process_sequence(seqwrt, &mut iter, &mut enzctxt, params, &mut unmap);
                    
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

pub fn generate_seqtable<R1: Read, R2: BufRead>(fasta: R1, tallymer: R2, params: &SeqTableParams, outfile: &str) {
    match params.mask {
        Some(_) => if params.strand_specific {
            generate_seqtable_ctxt(fasta, tallymer, params, EnzContextMaskedStrandSpecific::new(params), outfile)
        } else {
            generate_seqtable_ctxt(fasta, tallymer, params, EnzContextMasked::new(params), outfile)
        },
        None => generate_seqtable_ctxt(fasta, tallymer, params, EnzContextSimple::new(params.kmer_length, params.strand_specific), outfile),
    }
}

pub fn process_fasta(fasta_path: &str, tallymer_path: &OsStr, params: &SeqTableParams, outfile: &str) {
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