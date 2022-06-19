//!
//! Simple Enzyme Context - no masked bases
//!

use std::collections::VecDeque;
use super::{EnzContext, KmerIndex};

pub struct EnzContextSimple {
    cut_dna_value: u64,
    kmer_size: usize,
    revcomp_minus: bool,
    rev_dna_value: u64,
    mult: u64,

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
                if self.revcomp_minus {
                    self.rev_dna_value -= ( 3 - head as u64);
                    self.rev_dna_value /= 4;
                    self.mult /= 4;
                }
            }
        }
        self.buf.push_back(value);
        self.cut_dna_value *= 4;
        self.cut_dna_value += value as u64;
        if self.revcomp_minus {
            self.rev_dna_value += self.mult * (3 - value as u64);
            self.mult *= 4;
        }

        return self.buf.len() == self.kmer_size;
    }

    // public API
    pub fn new(kmer_size: u8, revcomp_minus: bool) -> EnzContextSimple {
        EnzContextSimple {
            cut_dna_value: 0u64,
            kmer_size: kmer_size as usize,
            revcomp_minus: revcomp_minus,
            rev_dna_value: 0u64,
            mult: 1,
            radix_power: 4u64.pow((kmer_size - 1) as u32),
            buf: VecDeque::new(),
        }
    }
}

impl EnzContext for EnzContextSimple {
    fn sequence_change(&mut self) {
        self.buf.clear();
        self.cut_dna_value = 0;
        self.mult = 1;
        self.rev_dna_value = 0;
    }

    fn add_base(&mut self, base: u8) -> KmerIndex {
        if self.update_cut_dna_value(base) {
            if self.revcomp_minus {
                KmerIndex { plus: Some(self.cut_dna_value as u32), minus: Some(self.rev_dna_value as u32) }
            } else {
                KmerIndex { plus: Some(self.cut_dna_value as u32), minus: Some(self.cut_dna_value as u32) }
            }
        } else {
            KmerIndex{ plus: None, minus: None }
        }
    }
}

#[cfg(test)]
mod tests {
    use fasta;
    use fasta::DNABases;
    use super::*;

    // AA 0   CA 4   GA 8   TA 12
    // AC     CC     GC     TC
    // AG     CG     GG     TG
    // AT     CT     GT     TT 15

    fn kmer_value(bases: Vec<DNABases>) -> Option<u32> {
        let mut res = 0;
        for base in bases {
            if base == DNABases::N {
                return None;
            }
            res *= 4;
            res += base as u32;
        }
        Some(res)
    }

    #[test]
    fn kmer_value_ok() {
        assert_eq!( Some(6), kmer_value(vec![DNABases::C, DNABases::G]));
        assert_eq!( Some(12), kmer_value(vec![DNABases::T, DNABases::A]));
        assert_eq!( None, kmer_value(vec![DNABases::C, DNABases::N]));
    }

    #[test]
    fn valid_first_index() {
        let mut simple = EnzContextSimple::new(2, false);
        assert_eq!( KmerIndex { plus: None, minus: None }, simple.add_base( fasta::DNABases::C as u8 ) );
        assert_eq!( KmerIndex {
            plus: kmer_value(vec![DNABases::C, DNABases::T]), // CT
            minus: kmer_value(vec![DNABases::C, DNABases::T]) // AG
        }, simple.add_base( fasta::DNABases::T as u8 ) );
    }

    #[test]
    fn valid_first_reverse_complement_index() {
        let mut simple = EnzContextSimple::new(2, true);
        assert_eq!( KmerIndex { plus: None, minus: None }, simple.add_base( fasta::DNABases::C as u8 ) );
        assert_eq!( KmerIndex {
            plus: kmer_value(vec![DNABases::C, DNABases::T]), // CT
            minus: kmer_value(vec![DNABases::A, DNABases::G]) // AG
        }, simple.add_base( fasta::DNABases::T as u8 ) );
    }

    #[test]
    fn valid_second_index() {
        let mut simple = EnzContextSimple::new(2, false);
        simple.add_base( fasta::DNABases::C as u8 );
        simple.add_base( fasta::DNABases::G as u8 );
        assert_eq!( KmerIndex {
            plus: kmer_value(vec![DNABases::G, DNABases::T]), // GT
            minus: kmer_value(vec![DNABases::G, DNABases::T]) // GT
        }, simple.add_base( fasta::DNABases::T as u8 ) );
    }

    #[test]
    fn valid_second_reverse_complement_index() {
        let mut simple = EnzContextSimple::new(2, true);
        simple.add_base( fasta::DNABases::C as u8 );
        simple.add_base( fasta::DNABases::G as u8 );
        assert_eq!( KmerIndex {
            plus: kmer_value(vec![DNABases::G, DNABases::T]), // GT
            minus: kmer_value(vec![DNABases::A, DNABases::C]) // AC
        }, simple.add_base( fasta::DNABases::T as u8 ) );

    }
}