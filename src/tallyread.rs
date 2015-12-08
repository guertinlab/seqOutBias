//!
//!	This modules is responsible for reading mappability information produced using tallymer.
//!
use std::io::prelude::*;
use std::ffi::OsString;
use std::io::Result;
use std::io::Error;
use std::io::ErrorKind;
use std::io::Lines;
use std::collections::BTreeMap;
use std::mem::swap;

/// Specifies if a position is unmappable in each strand.
#[derive(Copy, Clone)]
pub struct UnMapPosition {
    plus: bool,
    minus: bool,
}

pub struct UnMap<R: BufRead> {
    data: BTreeMap<u64, UnMapPosition>, // map from position to pair (unmappable in plus, unmappable in minus)
    data_next: BTreeMap<u64, UnMapPosition>, // will contain the data for the next sequence, if any 
    lines: Lines<R>,
    seqnumber: u64,
    next_seqnumber: u64,
}

fn parse_line(bytes: Vec<u8>) -> Result<(u64, u64, bool)> {
    let mut seq_idx = 0u64;
    let mut pos = 0u64;
    let mut idx = 0;
    
    while bytes[idx] >= b'0' && bytes[idx] <= b'9' {
        seq_idx = seq_idx * 10 + (bytes[idx] - b'0') as u64;
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
    while bytes[idx] >= b'0' && bytes[idx] <= b'9' {
        pos = pos * 10 + (bytes[idx] - b'0') as u64;
        idx += 1;
    }
    
    Ok((seq_idx, pos, is_minus))
}

impl<R: BufRead> UnMap<R> {
	
	/// Open mappability file, pre-loading data for the first sequence
	pub fn open(reader: R) -> Result<UnMap<R>> {
        let mut result = UnMap {
            data: BTreeMap::new(),
            data_next: BTreeMap::new(),
            lines: reader.lines(),
            seqnumber: -1,
            next_seqnumber: 0,
        };
        try!(result.read_next_sequence());
        
        return Ok(result);
	}
    
    fn insert_value(map: &mut BTreeMap<u64, UnMapPosition>, pos: u64, is_minus: bool) {
        if let Some(entry) = map.get_mut(&pos) {
            if is_minus {
                entry.minus = true;
            } else {
                entry.plus = true;
            }
            return;
        }
        
        // else (not written as an else to satisfy the borrow checker!)
        map.insert(pos, UnMapPosition { plus: !is_minus, minus: is_minus });
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
        while let Some(line) = self.lines.next() {
            let line = try!(line);
            let (seq, pos, is_minus) = try!(parse_line(line.into_bytes()));
            
            if seq == self.seqnumber {
                Self::insert_value(&mut self.data, pos, is_minus);
            } else {
                Self::insert_value(&mut self.data_next, pos, is_minus);
                self.next_seqnumber = seq;
                break;
            }
        }
        
        Ok(())
    }
    
    /// Query position
    ///
    /// Returns a pair of boolean values, which are true if the plus or minus strand read, respectively, is mappable at that position.
    pub fn is_mappable(&self, position: u64) -> UnMapPosition {
        match self.data.get(&position) {
            Some(&result) => result,
            None => UnMapPosition { plus: false, minus: false }, // if not in map, then it's not unmappable
        }
    }
}
