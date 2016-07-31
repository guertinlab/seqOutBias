//!
//!	Code to write a sequence table to disk.
//!
extern crate bincode;

use std::io::prelude::*;
use std::io::SeekFrom;
use std::io::Result;
use byteorder::{LittleEndian, WriteBytesExt};
use std::mem::size_of;
use bincode::rustc_serialize::{encode_into, encode, encoded_size};
use flate2::{Compression,Compress,Flush};
use std::cmp::max;

use ::seqtable::SeqTableParams;
use ::seqtable::SeqStore;

#[derive(Debug, RustcEncodable, RustcDecodable, PartialEq)]
pub struct SeqBlock {
    pub enc_size: u64,
    pub comp_size: u64,
    pub offset: u64,
}

#[derive(Debug, RustcEncodable, RustcDecodable, PartialEq)]
pub struct SeqInfo {
    pub name: String,
    pub length: u32,
    pub blocks: Vec<SeqBlock>,
}

pub struct SeqTableWriter<W: Write + Seek> {
    tailoffset: u64,
    writer: W,
    infotable: Vec<SeqInfo>,
    block_length: u32,
    max_buffer_size: u64,
    counts: Vec<(u64, u64, u64, u64)>,
}

impl<W: Write + Seek> SeqTableWriter<W> {
    pub fn new(mut writer: W, params: &SeqTableParams, block_length: u32) -> Result<SeqTableWriter<W>> {
        let blen = block_length;
        
        // write version number
        try!(writer.write_u8(super::TBL_VERSION));
        // write parameters to file
        try!(writer.write_u8(params.cut_length));
        // offsets are not required, but are stored anyway for reference
        try!(writer.write_u8(params.plus_offset));  
        try!(writer.write_u8(params.minus_offset));
        try!(writer.write_u16::<LittleEndian>(params.read_length));
        try!(writer.write_u32::<LittleEndian>(blen));
        
        // write temporary blank value to be filled in later
        // with offset of sequence table at end of file
        try!(writer.write_u64::<LittleEndian>(0));
        // same of max decoder size
        try!(writer.write_u64::<LittleEndian>(0));
        // same for counts table offset
        try!(writer.write_u64::<LittleEndian>(0));
        
        //
        let offset = size_of::<u32>() + 3 * size_of::<u64>() + 4 * size_of::<u8>() + size_of::<u16>();
        
        // allocate counts table
        let mut counts: Vec<(u64, u64, u64, u64)> = Vec::new();
        let nmer_count = 4u64.pow(params.cut_length as u32) + 1;
        
        for _ in 0..nmer_count {
            counts.push((0,0,0,0));
        }
        
        //
        Ok(SeqTableWriter{ 
            tailoffset: offset as u64, 
            writer: writer, 
            infotable: Vec::new(),
            block_length: blen,
            max_buffer_size: 0,
            counts: counts,
        })
    }
    
    pub fn create_sequence<'a>(&'a mut self, name: String) -> SequenceWriter<'a, W> {
        // create entry in info-table
        // return new sequence writter at current offset
        
        self.infotable.push( SeqInfo { name: name.clone(), length: 0, blocks: Vec::new() } );
        SequenceWriter {
            offset: &mut self.tailoffset,
            writer: &mut self.writer, 
            info: self.infotable.last_mut().unwrap(),
            block: Vec::new(),
            compressor: Compress::new(Compression::Best, false),
            output: vec![0u8; self.block_length as usize * size_of::<(u32,u32)>()],
            block_length: self.block_length,
            max_buffer_size: &mut self.max_buffer_size,
            counts: &mut self.counts,
        }
    }
}

impl<W: Write + Seek> Drop for SeqTableWriter<W> {
    fn drop(&mut self) {
        // at this point the writer is at the end of the file
        // and the tailoffset holds what will be the position
        // of the info table
        
        // write info table
        encode_into(&self.infotable, &mut self.writer, bincode::SizeLimit::Infinite).unwrap();
        
        // write counts table
        encode_into(&self.counts, &mut self.writer, bincode::SizeLimit::Infinite).unwrap();
        
        let counts_offset = self.tailoffset + encoded_size(&self.infotable);
        
        // seek to start & fill info table offset and max decoder size in header
        let offset = size_of::<u32>() + 4 * size_of::<u8>() + size_of::<u16>();
        self.writer.seek(SeekFrom::Start(offset as u64)).unwrap();
        self.writer.write_u64::<LittleEndian>(self.tailoffset as u64).unwrap();
        self.writer.write_u64::<LittleEndian>(self.max_buffer_size).unwrap();
        self.writer.write_u64::<LittleEndian>(counts_offset).unwrap();
    }
}

pub struct SequenceWriter<'a, W: 'a + Write> {
    offset: &'a mut u64,
    writer: &'a mut W,
    info: &'a mut SeqInfo,
    block: Vec<(u32, u32)>,
    compressor: Compress,
    output: Vec<u8>,
    block_length: u32,
    max_buffer_size: &'a mut u64,
    counts: &'a mut Vec<(u64, u64, u64, u64)>,
}

impl<'a, W: 'a + Write> SequenceWriter<'a, W> { 
    fn write_block(&mut self) -> Result<()> {
        /*for &(p, m) in &self.block {
            println!("write: ({},{})", p, m);
        }*/
        
        //   compress block and write to disk
        let binvec: Vec<u8> = encode(&self.block, bincode::SizeLimit::Infinite).unwrap();
        *self.max_buffer_size = max(*self.max_buffer_size, binvec.len() as u64);

        self.compressor.reset();
        let mut total = self.compressor.total_out();
        self.compressor.compress(&binvec, &mut self.output, Flush::Finish);
        total = self.compressor.total_out() - total;

        try!(self.writer.write(&self.output[0..(total as usize)]));
        
        //   add block to info-table && length
        self.info.blocks.push(SeqBlock { enc_size: binvec.len() as u64, comp_size: total, offset: *self.offset });
        *self.offset = *self.offset + total;
        
        self.info.length = self.info.length + self.block.len() as u32;
        self.block.clear();
        
        Ok(())
    }
}

impl<'a, W: 'a + Write> SeqStore for SequenceWriter<'a, W> {   
    fn write(&mut self, plus: u32, minus: u32) {
        // add to counts
        self.counts[plus as usize].0 += 1;
        self.counts[minus as usize].1 += 1;
        
        // add to block
        self.block.push((plus, minus));
        
        // store block if full
        if self.block.len() == self.block_length as usize {
                self.write_block().ok().expect("Failed to write block");
        }
    }
}

impl<'a, W: 'a + Write> Drop for SequenceWriter<'a, W> {
    fn drop(&mut self) {
        // write last block
        if self.block.len() > 0 {
            self.write_block().ok().expect("Failed to write block");
        }
    }
}
