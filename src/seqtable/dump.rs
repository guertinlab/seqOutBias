//!
//!	Code to write a sequence table to the terminal.
//!
use std::io::prelude::*;
use std::process::exit;
use std::fs::File;
use super::read::SeqTable;
use super::read::SequenceInfo;
use std::error::Error;

fn dump_header<R: Read + Seek>(table: &SeqTable<R>) {
    let params = table.params();
        
    println!("# kmer-size:    {}", params.kmer_length);
    println!("# plus-offset:  {}", params.plus_offset);
    println!("# minus-offset: {}", params.minus_offset);
    println!("# read-size:    {}", params.read_length);
    if let Some(ref mask) = params.mask {
        println!("# kmer-mask: {}", mask.iter().map(|&flag| if flag { 'N' } else { 'X' } ).collect::<String>());
    };
}

/// range can take the form: <chrom>:<start>-<end> or just <chrom>
fn decode_range(range: &str) -> (String, Option<(u32,u32)>) {
    let parts: Vec<&str> = range.split(':').collect();
    let chrom = parts[0].to_string();
    
    if parts.len() == 2 {
        let parts: Vec<&str> = parts[1].split('-').collect();
        
        if parts.len() == 2 {
            let start: u32 = parts[0].parse().expect("start coordinate");
            let end: u32 = parts[1].parse().expect("end coordinate");
            (chrom, Some((start, end)))
        } else {
            println!("Error: Malformed seqrange pattern {}.", range);
            exit(1);
        }
    } else {
        (chrom, None)
    }
}

/// Output only specified range
///
/// range can take the form: <chrom>:<start>-<end> or just <chrom>
pub fn dump_seqtable_range(filename: &str, range: &str) {
    // read
    let file = match File::open(filename) {
        Ok(value) => value,
        Err(err) => {
            println!("Error: failed to open sequence table file '{}': {}", filename, err.description());
            exit(1);
        }
    };
    let mut table = SeqTable::open(file).ok().expect("read store");
    
    // print header
    dump_header(&table);
    
    let seqinfos = table.sequences();

    // decode
    let coords = decode_range(range);
    
    for SequenceInfo { name, length } in seqinfos {
        if name.eq(&(coords.0)) {
            println!(">{} {}", name, length);
            let mut rdr = table.get_sequence(&name).ok().expect("read sequence");
            
            if let Some((start, end)) = coords.1 {
                for i in start..end {
                    let pair = rdr.get(i).unwrap();
                    println!("{}\t{}\t{}", i, pair.0, pair.1);
                }
            } else {
                for i in 0..length {
                    let pair = rdr.get(i).unwrap();
                    println!("{}\t{}\t{}", i, pair.0, pair.1);
                }
            }        
            
            break;
        }
    }
}

/// Output entire sequence table
pub fn dump_seqtable(filename: &str) {
    // read
    let file = match File::open(filename) {
        Ok(value) => value,
        Err(err) => {
            println!("Error: failed to open sequence table file '{}': {}", filename, err.description());
            exit(1);
        }
    };
    let mut table = SeqTable::open(file).ok().expect("read store");
    
    // print header
    dump_header(&table);
    
    let seqinfos = table.sequences();
    
    for SequenceInfo { name, length } in seqinfos {
        println!(">{} {}", name, length);
        
        let mut rdr = table.get_sequence(&name).ok().expect("read sequence");
        
        for i in 0..length {
            let pair = rdr.get(i).unwrap();
            println!("{}\t{}\t{}", i, pair.0, pair.1);
        }
    }
}
