//!
//!	This module and sub-modules contain the code to read and write the sequence cuts table using the mappability information
//!
use std::collections::VecDeque;
use tallyread::{UnMap,UnMapPosition};
use std::io::BufRead;
use std::cmp;

const TBL_VERSION : u8 = 2u8;

#[derive(Clone, Copy, Debug)]
pub struct SeqTableParams {
	pub cut_length: u8,
	pub plus_offset: u8,
	pub minus_offset: u8,
	pub read_length: u16,
}

pub trait SeqStore {
    fn write(&mut self, plus: u32, minus: u32);
}

mod write;
mod read;
mod dump;

// re-exports
pub use self::write::SeqTableWriter;
pub use self::write::SequenceWriter;
pub use self::read::SeqTable;
pub use self::read::SequenceInfo;
pub use self::dump::dump_seqtable;
pub use self::dump::dump_seqtable_range;

/// This buffer is used to translate between coordinate systems
/// Maps the n-mer table index values from the FASTA scan coordinates
/// to aligned read start coordinates.
pub struct SeqBuffer<'a, S: SeqStore, R: 'a + BufRead> {
    store: S,
    position: u32,
    written: u32, // TODO: consider storing this on SeqStore, i.e., add a method fn count(&self) -> u32
    plus_values: VecDeque<u32>,
    minus_values: VecDeque<u32>,
    plus_skip: u16,
    minus_skip: u16,
    common_skip: u16,
    unmap: &'a UnMap<R>,
}

impl<'a, S: SeqStore, R: BufRead> SeqBuffer<'a, S, R> {
    
    /// Create new sequence buffer which will store values into the supplied SeqStore instance
    pub fn new(mut store: S, params: SeqTableParams, unmap: &'a UnMap<R>) -> SeqBuffer<'a, S,R> {
        // Skipping
        //
        // If the aligned read position for each strand is before the start of the sequence, then
        // those values are not useful, i.e., they don't contribute to the output. So they should
        // be skipped, instead of stored in the value queues. 
        let plus_start = params.cut_length as i16 - 1 - params.plus_offset as i16;
        let minus_start = params.read_length as i16 - 1 + params.minus_offset as i16;
        let common_skip = if plus_start > 0 { // minus_start is always > 0 (assuming unsigned offsets)
          cmp::min(plus_start, minus_start)
        } else { 0 };
        
        // Padding
        //
        // Output needs padding when the first useful aligned read position for a strand is beyond the
        // start of the sequence. That can only happen for the plus strand (assuming unsigned offsets)
        // since it's the strand that can use the offset to get ahead of the cut-site region.
        let mut written = 0;
        if plus_start < 0 {
          let pad = -plus_start;
          for _ in 0..pad {
            store.write(0, 0);
          }
          written = pad as u32;
        }
        
        //
        SeqBuffer { 
          store: store,
          position: 0,
          written: written,
          plus_values: VecDeque::new(),
          minus_values: VecDeque::new(),
          plus_skip: if plus_start > 0 { plus_start as u16 } else { 0 },
          minus_skip: minus_start as u16,
          common_skip: common_skip as u16,
          unmap: unmap}
    }
    
    // Write values into underlying SeqStore, masking unmappable positions
    fn write(&mut self, plus_value: u32, minus_value: u32) {
      let UnMapPosition{ plus: unmap_plus, minus: unmap_minus } = self.unmap.is_unmappable(self.written);
            
      let idx_plus = if unmap_plus { 0 } else { plus_value };
      let idx_minus = if unmap_minus { 0 } else { minus_value };
      
      self.store.write(idx_plus, idx_minus);
      self.written += 1;
    }
    
    /// Push a new n-mer table index into buffer
    pub fn push(&mut self, table_index: u32) {
      if self.common_skip > 0 {
        self.common_skip -= 1;
        self.plus_skip -= 1;
        self.minus_skip -= 1;
      } else {
        if self.plus_skip == 0 && self.minus_skip == 0 {
          // obtain plus strand value at current position
          let plus_value = if self.plus_values.len() > 0 {
            self.plus_values.push_front(table_index);
            self.plus_values.pop_back().unwrap()
          } else { table_index };
          
          // obtain minus strand value at current position
          let minus_value = if self.minus_values.len() > 0 {
            self.minus_values.push_front(table_index);
            self.minus_values.pop_back().unwrap()
          } else { table_index };
          
          // write pairs to store
          self.write(plus_value, minus_value);
        } else {
          // Store values until the remaining strand has caught up
          if self.minus_skip > 0 {
            self.plus_values.push_front(table_index);
            self.minus_skip -= 1;
          }
          if self.plus_skip > 0 {
            self.minus_values.push_front(table_index);
            self.plus_skip -= 1;
          }
        }
      }
      
      // update sequence position
      self.position += 1;
    }
}

impl<'a, S: SeqStore, R: BufRead> Drop for SeqBuffer<'a, S, R> {
  
  // On Drop, complete the output as needed using remaining buffer content
  fn drop(&mut self) {
    let seq_length = self.position;
    
    // Output needs padding up to the length of the sequence
    while self.written < seq_length {
      let plus_value = self.plus_values.pop_back().unwrap_or(0);
      let minus_value = self.minus_values.pop_back().unwrap_or(0);
      
      self.write(plus_value, minus_value);
    }
  }
}
