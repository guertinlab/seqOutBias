//!
//! Masked Enzyme Context - kmer can have masked positions
//!

use std::collections::VecDeque;
use super::{EnzContext, KmerIndex};
use seqtable::SeqTableParams;

// Masked version
pub struct EnzContextMasked {
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
