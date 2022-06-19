//!
//! Masked Enzyme Context with strand specific behaviour - kmer can have masked positions
//!

use std::collections::VecDeque;
use super::{EnzContext, KmerIndex};
use seqtable::SeqTableParams;

// Strand Specific Masked version
pub struct EnzContextMaskedStrandSpecific {
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
