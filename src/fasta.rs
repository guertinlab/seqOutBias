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

/// Observed kmer on the plus and minus strands
/// Unless the strand-specific flag is active, both values are equal
/// and match the plus strand readout.
struct KmerIndex {
    plus: Option<u32>,
    minus: Option<u32>
}

/// Holds information on the current enzyme cut region
trait EnzContext {
    fn sequence_change(&mut self);
    fn add_base(&mut self, base: u8) -> KmerIndex;
}

struct EnzContextSimple {
    cut_dna_value: u64,
    kmer_size: usize,
    revcomp_minus: bool,
    
    // auxiliary data
    radix_power: u64,
    buf: VecDeque<u8>,
}

impl EnzContextSimple {
    // private API
    fn update_cut_dna_value(&mut self, value: u8) -> bool {
        if value == 4 {
            self.buf.clear();
            self.cut_dna_value = 0;
            return false;
        }

        if self.buf.len() == self.kmer_size {
            if let Some(head) = self.buf.pop_front() {
                self.cut_dna_value -= self.radix_power * (head as u64);
            }
        }
        self.buf.push_back(value);
        self.cut_dna_value *= 4;
        self.cut_dna_value += value as u64;
    
        return self.buf.len() == self.kmer_size;
    }
        
    // public API
    pub fn new(kmer_size: u8, revcomp_minus: bool) -> EnzContextSimple {
        EnzContextSimple {
            cut_dna_value: 0u64,
            kmer_size: kmer_size as usize,
            revcomp_minus: revcomp_minus,
            radix_power: 4u64.pow((kmer_size - 1) as u32),
            buf: VecDeque::new(),
        }
    }
}

impl EnzContext for EnzContextSimple {
    fn sequence_change(&mut self) {
        self.buf.clear();
        self.cut_dna_value = 0;
    }
    
    fn add_base(&mut self, base: u8) -> KmerIndex {
        if self.update_cut_dna_value(base) {
            if self.revcomp_minus {
                let revcomp = 4u32.pow( self.kmer_size as u32 ) - self.cut_dna_value as u32 - 1;
                KmerIndex { plus: Some(self.cut_dna_value as u32), minus: Some(revcomp) }
            } else {
                KmerIndex { plus: Some(self.cut_dna_value as u32), minus: Some(self.cut_dna_value as u32) }
            }
        } else {
            KmerIndex{ plus: None, minus: None }
        }
    }
}

// Masked version
struct EnzContextMasked {
    kmer_length: usize,
    unmasked_count: u32,
    buf: VecDeque<u8>,
    mask: Vec<bool>,
}

impl EnzContextMasked {
    pub fn new(params: &SeqTableParams) -> EnzContextMasked {
        EnzContextMasked {
            kmer_length: params.kmer_length as usize,
            unmasked_count: params.unmasked_count as u32,
            buf: VecDeque::with_capacity(params.unmasked_count as usize),
            mask: params.mask.as_ref().unwrap().clone(),
        }
    }

    fn nmer_index(&self) -> KmerIndex {
        let mut idx = 0;
        let mut mult = 4u64.pow(self.unmasked_count - 1);

        for (base, flag) in self.buf.iter().zip(&self.mask) {
            if *flag {
                if *base == 4 { return KmerIndex{ plus: None, minus: None }; } // N in a non-masked position

                idx = idx + *base as u64 * mult;
                mult = mult / 4;
            }
        }
        KmerIndex{ plus: Some( idx as u32 ), minus: Some( idx as u32 ) }
    }
}

impl EnzContext for EnzContextMasked {
    fn sequence_change(&mut self) {
        self.buf.clear();
    }
    
    fn add_base(&mut self, base: u8) -> KmerIndex {
        if self.buf.len() == self.kmer_length { self.buf.pop_front(); }
        self.buf.push_back(base);
        if self.buf.len() == self.kmer_length {
            self.nmer_index()
        } else {
            KmerIndex{ plus: None, minus: None }
        }
    }
}

// Strand Specific Masked version
struct EnzContextMaskedStrandSpecific {
    kmer_length: usize,
    unmasked_count: u32,
    buf: VecDeque<u8>,
    mask_plus: Vec<bool>,
    mask_minus: Vec<bool>,
}

impl EnzContextMaskedStrandSpecific {
    pub fn new(params: &SeqTableParams) -> EnzContextMaskedStrandSpecific {
        let mut rev_mask = params.mask.as_ref().unwrap().clone();
        rev_mask.reverse();
        EnzContextMaskedStrandSpecific {
            kmer_length: params.kmer_length as usize,
            unmasked_count: params.unmasked_count as u32,
            buf: VecDeque::with_capacity(params.unmasked_count as usize),
            mask_plus: params.mask.as_ref().unwrap().clone(),
            mask_minus: rev_mask,
        }
    }

    fn nmer_index(&self) -> KmerIndex {
        let mut idx_plus = 0;
        let mut idx_minus = 0;
        let mut plus_is_some = true;
        let mut minus_is_some = true;
        let mut mult = 4u64.pow(self.unmasked_count - 1);

        // compute plus kmer
        for (base, flag) in self.buf.iter().zip(&self.mask_plus) {
            if *flag {
                if *base == 4 { plus_is_some = false; break } // N in a non-masked position

                idx_plus = idx_plus + *base as u64 * mult;
                mult = mult / 4;
            }
        }
        // compute minus kmer - reverse complement
        mult = 1;
        for (base, flag) in self.buf.iter().zip(&self.mask_minus) {
            if *flag {
                if *base == 4 { minus_is_some = false; break } // N in a non-masked position

                idx_minus = idx_minus + ( 3 - *base ) as u64 * mult;
                mult = mult * 4;
            }
        }
        KmerIndex {
            plus: if plus_is_some { Some(idx_plus as u32) } else { None },
            minus: if minus_is_some { Some(idx_minus as u32) } else { None }
        }
    }
}

impl EnzContext for EnzContextMaskedStrandSpecific {
    fn sequence_change(&mut self) {
        self.buf.clear();
    }

    fn add_base(&mut self, base: u8) -> KmerIndex {
        if self.buf.len() == self.kmer_length { self.buf.pop_front(); }
        self.buf.push_back(base);
        if self.buf.len() == self.kmer_length {
            self.nmer_index()
        } else {
            KmerIndex{ plus: None, minus: None }
        }
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
            b'a' | b'A' => { store_base!(enzctxt, buf, 0); seqpos += 1; },
            b'c' | b'C' => { store_base!(enzctxt, buf, 1); seqpos += 1; },
            b'g' | b'G' => { store_base!(enzctxt, buf, 2); seqpos += 1; },
            b't' | b'T' => { store_base!(enzctxt, buf, 3); seqpos += 1; },
            b'n' | b'N' => { store_base!(enzctxt, buf, 4); seqpos += 1; },
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