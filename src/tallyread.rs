//!
//!	This module is responsible for reading mappability information produced using tallymer.
//!
use std::io::prelude::*;
use std::io::Result;
use std::io::Error;
use std::io::ErrorKind;
use std::mem::swap;

/// Specifies if a position is unmappable in each strand.
#[derive(Copy, Clone)]
pub struct UnMapPosition {
    pub plus: bool,
    pub minus: bool,
}

struct UnMapCell {
    position: u32,
    data: UnMapPosition,
}

pub struct UnMap<R: BufRead> {
    data: Vec<UnMapCell>, // map from position to pair (unmappable in plus, unmappable in minus)
    data_next: Vec<UnMapCell>, // will contain the data for the next sequence, if any 
    reader: R,
    line: Vec<u8>,
    seqnumber: u32,
    next_seqnumber: u32,
}

fn parse_line(bytes: &[u8]) -> Result<(u32, u32, bool)> {
    let mut seq_idx = 0u32;
    let mut pos = 0u32;
    let mut idx = 0;
        
    while bytes[idx] >= b'0' && bytes[idx] <= b'9' {
        seq_idx = seq_idx * 10 + (bytes[idx] - b'0') as u32;
        idx += 1;
    }
    
    if bytes[idx] != b'\t' {
        return Err(Error::new(ErrorKind::Other, "unexpected byte in tallymer line"));
    }
    
    // parse sign
    idx += 1;
    let is_minus = match bytes[idx] {
        b'-' => true,
        b'+' => false,
        _ => return Err(Error::new(ErrorKind::Other, "unexpected byte in tallymer line"))
    };
    idx += 1;
    
    // parse position
    while idx < bytes.len() && bytes[idx] >= b'0' && bytes[idx] <= b'9' {
        pos = pos * 10 + (bytes[idx] - b'0') as u32;
        idx += 1;
    }
    
    Ok((seq_idx, pos, is_minus))
}

impl<R: BufRead> UnMap<R> {
	
	/// Open mappability file, pre-loading data for the first sequence
	pub fn open(reader: R) -> Result<UnMap<R>> {
        let mut result = UnMap {
            data: Vec::new(),
            data_next: Vec::new(),
            reader: reader,
            line: Vec::new(),
            seqnumber: 0,
            next_seqnumber: 0,
        };
        try!(result.read_sequence_data());
        
        return Ok(result);
	}
    
    fn insert_value(map: &mut Vec<UnMapCell>, pos: u32, is_minus: bool) {
        match map.last_mut() {
            Some(ref mut entry) if entry.position == pos => {
                if is_minus {
                    entry.data.minus = true;
                } else {
                    entry.data.plus = true;
                }

                return;
            },
            _ => {},
        }
        map.push(UnMapCell{
                position: pos,
                data: UnMapPosition {
                    plus: !is_minus,
                    minus: is_minus,
                }
        });
    }
	
    fn read_sequence_data(&mut self) -> Result<()> {
        // read until either no more data or if find a row for the next sequence
        self.line.clear();
        while self.reader.read_until(b'\n', &mut self.line).unwrap() > 0 {
            let (seq, pos, is_minus) = try!(parse_line(&self.line));
            
            if seq == self.seqnumber {
                Self::insert_value(&mut self.data, pos, is_minus);
            } else {
                Self::insert_value(&mut self.data_next, pos, is_minus);
                self.next_seqnumber = seq;
                break;
            }

            // clear buffer
            self.line.clear();
        }
        Ok(())
    }
    
	/// Read the information for the next sequence from disk
	pub fn read_next_sequence(&mut self) -> Result<()> {
        // if we haven't reached the next sequence number yet
		if self.next_seqnumber > self.seqnumber + 1 {
            self.seqnumber = self.seqnumber + 1;
            self.data.clear();
            return Ok(());
        } else if self.seqnumber == self.next_seqnumber {
            // there is no more data to read
            self.data.clear();
            return Ok(());
        }
        self.seqnumber = self.next_seqnumber;
        
        // prepare
        swap(&mut self.data, &mut self.data_next);
        self.data_next.clear();
        
        // read until either no more data or if find a row for the next sequence
        try!(self.read_sequence_data());
        
        Ok(())
    }
    
    /// Query position
    ///
    /// Returns a pair of boolean values, which are true if the plus or minus strand read, respectively, is mappable at that position.
    pub fn is_unmappable(&self, position: u32) -> UnMapPosition {
        match self.data.binary_search_by( |probe| probe.position.cmp(&position) ) {
            Ok(index) => {
                unsafe {
                    self.data.get_unchecked(index).data
                }
            },
            _ => UnMapPosition { plus: false, minus: false }, // if not in map, then it's not unmappable
        }
    }

    /// get current sequence number
    pub fn sequence_number(&self) -> u32 {
        self.seqnumber
    }
}
