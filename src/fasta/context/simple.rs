//!
//! Simple Enzyme Context - no masked bases
//!

use std::collections::VecDeque;
use super::{EnzContext, KmerIndex};

pub struct EnzContextSimple {
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

