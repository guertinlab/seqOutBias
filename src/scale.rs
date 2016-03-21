//!
//! Code to create a scaled read-count track 
//!
use htslib::bam;
use htslib::bam::Read;
use htslib::bam::Reader;
use htslib::bam::Records;
use std::fs::File;
use std::io::Read as ioRead;
use std::io::Write as ioWrite;
use std::io::Seek;
use std::io::Error;
use std::process::exit;
use std::iter::Peekable;
use seqtable::{SeqTable,SequenceInfo};
use std::collections::BTreeMap;

#[derive(Debug)]
pub struct PileUp {
    chroms: Vec<String>,
    counts: Vec<BTreeMap<u32, (f64, f64)>>,
}

impl PileUp {
    
    fn new(sinfos: Vec<SequenceInfo>) -> PileUp {
        let mut chroms = Vec::new();
        let mut counts = Vec::new();
        
        for sinfo in sinfos {
            chroms.push(sinfo.name.clone());
            counts.push(BTreeMap::new());
        }
        
        PileUp { chroms: chroms, counts: counts }
    }
    
    fn add_data<R: ioRead+Seek>(&mut self, table: &mut SeqTable<R>, bamrecs: &mut Peekable<Records<Reader>>, tid: &mut i32, map: &Vec<usize>, rlen: usize, minqual: u8, counts: &Vec<(u64, u64, u64, u64)>) -> bool {
        let sidx = map[*tid as usize];
        let mut rdr = table.get_sequence_by_idx(sidx).ok().expect("read sequence");
    
        loop {
            // check if we changed sequence
            match bamrecs.peek() {
                Some(&Ok(ref rec)) => if rec.tid() != *tid { *tid = rec.tid(); return true; },
                Some(&Err(_)) => return false,
                None => return false,
            }
            
            // if not count position
            if let Some(Ok(record)) = bamrecs.next() {
                if !record.is_unmapped() && record.seq().len() == rlen && record.mapq() >= minqual {
                    let (plus_idx, minus_idx) = rdr.get(record.pos() as u32).unwrap();
                    
                    if record.is_reverse() {
                        if minus_idx == 0 {
                            /* no data */
                        } else {
                            let s_minus = counts[minus_idx as usize].1 as f64;
                            let b_minus = counts[minus_idx as usize].3 as f64;
                            if s_minus > 0f64 {
                                let inc = b_minus / s_minus;
                                self.counts[sidx as usize].entry(record.pos() as u32).or_insert((0f64, 0f64)).1 += inc;
                            } /* else no data */
                        }
                    } else {
                        if plus_idx == 0 {
                            /* no data */
                        } else {
                            let s_plus = counts[plus_idx as usize].0 as f64;
                            let b_plus = counts[plus_idx as usize].2 as f64;
                            if s_plus > 0f64 {
                                let inc = b_plus / s_plus;
                                self.counts[sidx as usize].entry(record.pos() as u32).or_insert((0f64, 0f64)).0 += inc;
                            } /* else no data */
                        }
                    } 
                }
            } 
        }
    }
    
    pub fn write_bed(&self, filename: &str, stranded: bool) -> Result<(),Error> {
        // open new file
        let mut f = try!(File::create(filename));
        
        // write data
        for i in 0..self.chroms.len() {
            let chrom : &str = &self.chroms[i];
            
            for (pos, value) in self.counts[i].iter() {
                
                if stranded {
                    try!(write!(f, "{}\t{}\t{}\t{}\n", chrom, pos, pos + 1, value.0));
                    try!(write!(f, "{}\t{}\t{}\t{}\n", chrom, pos, pos + 1, -value.1));
                } else {
                    try!(write!(f, "{}\t{}\t{}\t{}\n", chrom, pos, pos + 1, value.0 + value.1));
                }
            }
        }
        Ok(())
    }
}

pub fn scale(seqfile: &str, counts: Vec<(u64, u64, u64, u64)>, bamfile: String, minqual: u8) -> PileUp {
    // read
    let file = File::open(seqfile).ok().expect("read file");
    let mut table = match SeqTable::open(file) {
        Ok(value) => value,
        Err(e) => {
            println!("Error:tabulate: {}", e.to_string()); 
            exit(1); 
        },
    };
    
    let rlen = table.params.read_length as usize;
    let seqinfos = table.sequences();

    let bam = bam::Reader::new(&bamfile).ok().expect("Error opening bam.");
    let names = bam.header.target_names();
    let mut cur_tid = 0;
    
    // map BAM tid's to SeqTable idx's
    let map : Vec<usize> = names.iter().map(|&id| {
        let chrom = String::from_utf8_lossy(id);
        for i in 0..seqinfos.len() {
            if seqinfos[i].name.eq(&chrom) {
                return i;
            }
        }
        println!("Error: Unknown sequence name in BAM: {}", chrom);
        exit(1);
    }
    ).collect();
    
    // reads
    let mut iter = bam.records().peekable();
    let mut pileup = PileUp::new(seqinfos);
    
    while pileup.add_data(&mut table, &mut iter, &mut cur_tid, &map, rlen, minqual, &counts) {}
    
    pileup
}
