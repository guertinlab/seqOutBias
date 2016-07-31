//!
//!	Code to read a sequence table from disk.
//!
extern crate bincode;

use std::io::prelude::*;
use std::io::SeekFrom;
use std::io::Result;
use std::io::Error;
use std::io::ErrorKind;
use byteorder::{LittleEndian, ReadBytesExt};
use bincode::rustc_serialize::{decode_from, decode};
use flate2::{Decompress,Flush};

use super::SeqTableParams;
use super::write::SeqInfo;
use super::write::SeqBlock;

#[derive(Debug)]
pub struct SeqTable<R: Read + Seek> {
    pub params: SeqTableParams,
    block_length: u32,
    infotable: Vec<SeqInfo>,
    reader: R,
    dec_buffer: Vec<u8>,
    read_buffer: Vec<u8>,
    counts_offset: u64,
}

pub struct SequenceInfo {
    pub name: String,
    pub length: u32,
}

impl<R: Read + Seek> SeqTable<R> {
    pub fn open(mut reader: R) -> Result<SeqTable<R>> {
        try!(reader.seek(SeekFrom::Start(0)));
        // load version
        let version = try!(reader.read_u8());
        if version != super::TBL_VERSION {
            return Err(Error::new(ErrorKind::InvalidData, format!("Incompatible file version {}, expected {}.", version, super::TBL_VERSION)));
        }
        // load parameters
        let mut params = SeqTableParams {
            cut_length: try!(reader.read_u8()),
            plus_offset: try!(reader.read_u8()),
            minus_offset: try!(reader.read_u8()),
            read_length: try!(reader.read_u16::<LittleEndian>()),
            mask: None,
            unmasked_count: 0,
        };
        params.unmasked_count = params.cut_length;
         
        // load block size
        let blen = try!(reader.read_u32::<LittleEndian>());
        // load info table offset
        let offset = try!(reader.read_u64::<LittleEndian>());
        // load buffer size
        let bufsize = try!(reader.read_u64::<LittleEndian>()) as usize;
        // load counts table offset
        let counts_offset = try!(reader.read_u64::<LittleEndian>());
        
        //println!("blen: {}, offset: {}, bufsize: {}", blen, offset, bufsize);
        
        // load info table
        try!(reader.seek(SeekFrom::Start(offset)));
        let infotable = (decode_from(&mut reader, bincode::SizeLimit::Infinite)).unwrap(); // TODO: fix
        
        Ok(SeqTable {
            params: params,
            block_length: blen,
            infotable: infotable,
            reader: reader,
            dec_buffer: vec![0u8; bufsize],
            read_buffer: vec![0u8; bufsize],
            counts_offset: counts_offset,
        })
    }
    
    pub fn get_sequence<'a>(&'a mut self, name: &str) -> Result<SeqReader<'a, R>> {
        // locate sequence info
        match self.infotable.iter().find(|&x| x.name == name) {
            Some(iref) => {
                match iref.blocks.get(0) {
                    Some(block) => {
                        try!(self.reader.seek(SeekFrom::Start(block.offset)));
                        Ok(SeqReader { 
                            reader: &mut self.reader,
                            block_length: self.block_length,
                            info: &iref,
                            block: None,
                            block_idx: 0,
                            dec_buffer: &mut self.dec_buffer,
                            read_buffer: &mut self.read_buffer,
                        })
                    },
                    None => Err(Error::new(ErrorKind::NotFound, "sequence not found")), 
                }
                
            },
            None => Err(Error::new(ErrorKind::NotFound, "sequence not found")),
        }
    }
    
    pub fn get_sequence_by_idx<'a>(&'a mut self, idx: usize) -> Result<SeqReader<'a, R>> {
        if self.infotable.len() <= idx {
            Err(Error::new(ErrorKind::NotFound, "sequence not found"))
        } else {
            match self.infotable[idx].blocks.get(0) {
                Some(block) => {
                    try!(self.reader.seek(SeekFrom::Start(block.offset)));
                    Ok(SeqReader { 
                        reader: &mut self.reader,
                        block_length: self.block_length,
                        info: &self.infotable[idx],
                        block: None,
                        block_idx: 0,
                        dec_buffer: &mut self.dec_buffer,
                        read_buffer: &mut self.read_buffer,
                    })
                },
                None => Err(Error::new(ErrorKind::NotFound, "sequence not found")), 
            }
        }
    }
    
    pub fn params(&self) -> &SeqTableParams {
        &self.params
    }
    
    pub fn equivalent(&self, _: &str, params: &SeqTableParams) -> bool {
        self.params.read_length == params.read_length &&
        self.params.cut_length == params.cut_length &&
        self.params.plus_offset == params.plus_offset &&
        self.params.minus_offset == params.minus_offset
    }
    
    pub fn sequences(&self) -> Vec<SequenceInfo> {
        self.infotable.iter().map(|info| 
            SequenceInfo{ name: info.name.clone(), length: info.length }
        ).collect()
    }
    
    pub fn counts(&mut self) -> Result<Vec<(u64, u64, u64, u64)>> {
        try!(self.reader.seek(SeekFrom::Start(self.counts_offset)));
        let counts: Vec<(u64, u64, u64, u64)> = (decode_from(&mut self.reader, bincode::SizeLimit::Infinite)).unwrap();
        Ok(counts)
    }
    
    pub fn len(&self) -> usize {
        self.infotable.len()
    }
}

#[derive(Debug)]
pub struct SeqReader<'a, R: 'a + Read + Seek> {
    reader: &'a mut R,
    block_length: u32,
    info: &'a SeqInfo,
    block: Option<Vec<(u32, u32)>>,
    block_idx: usize,
    dec_buffer: &'a mut Vec<u8>,
    read_buffer: &'a mut Vec<u8>,
}

impl<'a, R: 'a + Read + Seek> SeqReader<'a, R> {
    
    fn read_block(&mut self, block_info: &SeqBlock) -> Result<Vec<(u32,u32)>> {
        
        // read from disk (compressed)
        let buf = &mut self.read_buffer[0..block_info.comp_size as usize];
        let mut decompressor = Decompress::new(false);
        
        try!(self.reader.seek(SeekFrom::Start(block_info.offset)));
        try!(self.reader.read(buf));
        
        // uncompress
        let bufout = &mut self.dec_buffer[0..block_info.enc_size as usize];
        decompressor.decompress(&buf, bufout, Flush::Finish).ok().expect("decompressing block");

        // decode back into vector
        let block: Vec<(u32, u32)> = decode(bufout).ok().expect("decode block");
        
        return Ok(block);
    }
    
    pub fn vir_get(&mut self, position: i32) -> Result<(u32, u32)> { 
        if position < 0 { Ok((0, 0)) }
        else { self.get(position as u32) }
    }

    pub fn get(&mut self, position: u32) -> Result<(u32, u32)> {
        let idx = (position / self.block_length) as usize;
        
        //println!("get: pos: {} idx: {} blen: {} n_blocks: {} length: {}", position, idx, self.block_length, self.info.blocks.len(), self.info.length);
        
        let block = if self.block.is_some() && self.block_idx == idx {
            self.block.as_ref().unwrap()
        } else {
            match self.info.blocks.get(idx) {
                Some(block_info) => {
                    let aux = try!(self.read_block(block_info));
                    self.block = Some(aux);
                    self.block_idx = idx;
                    self.block.as_ref().unwrap()
                },
                None => return Err(Error::new(ErrorKind::AddrNotAvailable, "index out of range")),
            }
        };
        
        // get value
        let idx2 = position as usize - idx * self.block_length as usize; 
        Ok(block[idx2])
    }
}
