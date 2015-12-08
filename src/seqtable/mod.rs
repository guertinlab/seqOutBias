//!
//!	This modules is contains the code to read and write the sequence cuts table using the mappability information
//!
use std::collections::VecDeque;

#[derive(Clone, Copy, Debug)]
pub struct SeqTableParams {
	pub cut_length: u8,
	pub plus_offset: u8, 
	pub minus_offset: u8,
	pub read_length: u8,
}

pub trait SeqStore {
    fn write(&mut self, plus: u16, minus: u16);
}

mod write;
mod read;

// re-exports
pub use self::write::SeqTableWriter;
pub use self::read::SeqTable;

/// This buffer is used to translate between coordinate systems
/// Maps the n-mer table index values from the FASTA scan coordinates
/// to aligned read start coordinates.
pub struct SeqBuffer<S: SeqStore> {
    store: S,
    length: usize,
    plus_values: VecDeque<u16>,
}

impl<S: SeqStore> SeqBuffer<S> {
    
    /// Create new sequence buffer which will store values into the supplied SeqStore instance
    pub fn new(store: S, params: SeqTableParams) -> SeqBuffer<S> {
        let buf_len = params.read_length + params.minus_offset - params.cut_length + params.plus_offset;
        let mut queue = VecDeque::new();
        
        // pad buffer
        // for the first OP positions there are no values to insert because the cut-site falls off the edge of the sequence
        if params.plus_offset > 0 {
            for _ in 0..params.plus_offset {
                queue.push_front(0);
            }
        }
        
        SeqBuffer { store: store, length: buf_len as usize, plus_values: queue }
    }
    
    /// Push a new n-mer table index into buffer
    pub fn push(&mut self, table_index: u16) {
        self.plus_values.push_front(table_index);
        
        if self.plus_values.len() > self.length {
            let value = self.plus_values.pop_back().unwrap();
            
            self.store.write(value, table_index);
        }
    }
}