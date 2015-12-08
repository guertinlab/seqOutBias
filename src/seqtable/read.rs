//!
//!	Code to write a sequence table to disk.
//!
extern crate bincode;

use std::io::prelude::*;
use std::io::SeekFrom;
use std::io::Result;
use std::io::Error;
use std::io::ErrorKind;
use byteorder::{ByteOrder, LittleEndian, ReadBytesExt};
use bincode::rustc_serialize::{decode_from, decode};
use flate2::{Decompress,Flush};

use super::SeqTableParams;
use super::write::SeqInfo;
use super::write::SeqBlock;

#[derive(Debug)]
pub struct SeqTable<R: Read + Seek> {
    params: SeqTableParams,
    block_length: u64,
    infotable: Vec<SeqInfo>,
    reader: R,
    dec_buffer: Vec<u8>,
    read_buffer: Vec<u8>,
}

impl<R: Read + Seek> SeqTable<R> {
    pub fn open(mut reader: R) -> Result<SeqTable<R>> {
        try!(reader.seek(SeekFrom::Start(0)));
        // load parameters
        let params = SeqTableParams {
            cut_length: try!(reader.read_u8()),
            plus_offset: try!(reader.read_u8()),
            minus_offset: try!(reader.read_u8()),
            read_length: try!(reader.read_u8()),
        };
        // load block size
        let blen = try!(reader.read_u64::<LittleEndian>());
        // load info table offset
        let offset = try!(reader.read_u64::<LittleEndian>());
        // load buffer size
        let bufsize = try!(reader.read_u64::<LittleEndian>()) as usize;
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
}

#[derive(Debug)]
struct SeqReader<'a, R: 'a + Read + Seek> {
    reader: &'a mut R,
    block_length: u64,
    info: &'a SeqInfo,
    block: Option<Vec<(u16, u16)>>,
    block_idx: usize,
    dec_buffer: &'a mut Vec<u8>,
    read_buffer: &'a mut Vec<u8>,
}

impl<'a, R: 'a + Read + Seek> SeqReader<'a, R> {
    
    fn read_block(&mut self, block_info: &SeqBlock) -> Result<Vec<(u16,u16)>> {
        
        // read from disk (compressed)
        let buf = &mut self.read_buffer[0..block_info.comp_size as usize];
        let mut decompressor = Decompress::new(false);
        
        try!(self.reader.seek(SeekFrom::Start(block_info.offset)));
        try!(self.reader.read(buf));
        
        // uncompress
        let bufout = &mut self.dec_buffer[0..block_info.enc_size as usize];
        decompressor.decompress(&buf, bufout, Flush::Finish).ok().expect("decompressing block");

        // decode back into vector
        let block: Vec<(u16, u16)> = decode(bufout).ok().expect("decode block");
        
        return Ok(block);
    }
    
    pub fn get(&mut self, position: u64) -> Result<(u16, u16)> {
        let idx = (position / self.block_length) as usize;
        
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
