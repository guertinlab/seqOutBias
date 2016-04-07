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
    minus_shift: i32,
    no_scale: bool,
}

impl PileUp {
    
    fn new(sinfos: Vec<SequenceInfo>, minus_shift: i32, no_scale: bool) -> PileUp {
        let mut chroms = Vec::new();
        let mut counts = Vec::new();
        
        for sinfo in sinfos {
            chroms.push(sinfo.name.clone());
            counts.push(BTreeMap::new());
        }
        
        PileUp { chroms: chroms, counts: counts, minus_shift: minus_shift, no_scale: no_scale}
    }
    
    fn add_data<R: ioRead+Seek>(&mut self, table: &mut SeqTable<R>, bamrecs: &mut Peekable<Records<Reader>>, tid: &mut i32, map: &Vec<usize>, minqual: u8, scale: &Vec<(f64, f64)>) -> bool {
        // skip unmapped sequences (tid = -1)
        if *tid < 0 {
            match bamrecs.next() {
                Some(Ok(ref rec)) => { *tid = rec.tid(); return true; },
                Some(Err(_)) => return false,
                None => return false,
            }
        }
        
        let sidx = map[*tid as usize];
        let rlen = table.params.read_length as usize;
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
                            let inc = if self.no_scale { 1f64 } else { scale[minus_idx as usize].1 };
                            let minus_pos = ((record.pos() as u32 + rlen as u32 - 1u32) as i32 + self.minus_shift) as u32;
                            self.counts[sidx as usize].entry(minus_pos).or_insert((0f64, 0f64)).1 += inc;
                        }
                    } else {
                        if plus_idx == 0 {
                            /* no data */
                        } else {
                            let inc = if self.no_scale { 1f64 } else { scale[plus_idx as usize].0 };
                            self.counts[sidx as usize].entry(record.pos() as u32).or_insert((0f64, 0f64)).0 += inc;
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
                    if value.0 > 0f64 {
                        try!(write!(f, "{}\t{}\t{}\t.\t{}\t+\n", chrom, pos, pos + 1, value.0));
                    }
                    if value.1 > 0f64 {
                        try!(write!(f, "{}\t{}\t{}\t.\t{}\t-\n", chrom, pos, pos + 1, -value.1));
                    }
                } else {
                    try!(write!(f, "{}\t{}\t{}\t.\t{}\n", chrom, pos, pos + 1, value.0 + value.1));
                }
            }
        }
        Ok(())
    }
}

fn scale_factor(exp: u64, etotal: u64, obs: u64, ototal: u64) -> f64 {
    let fexp = if etotal > 0 { exp as f64 / etotal as f64 } else { 0f64 };
    let fobs = if ototal > 0 { obs as f64 / ototal as f64 } else { 0f64 };
    if fobs > 0f64 { fexp / fobs } else { 0f64 }
}

fn compute_scale_factors(counts: &Vec<(u64, u64, u64, u64)>) -> Vec<(f64, f64)> {
    // compute totals
    let totals = counts.iter().fold((0u64, 0u64, 0u64, 0u64),
        |acc, &(sp, sm, bp, bm)| ( acc.0 + sp, acc.1 + sm, acc.2 + bp, acc.3 + bm));
    
    // compute scale
    // x = Obs * ExpFreq / ObsFreq
    // scale = ExpFreq / ObsFreq
    counts.iter().map(|&(sp, sm, bp, bm)| {
        ( scale_factor(sp, totals.0, bp, totals.2),
          scale_factor(sm, totals.1, bm, totals.3)
        )
    }).collect()
}

pub fn scale(seqfile: &str, counts: Vec<(u64, u64, u64, u64)>, bamfile: String, minqual: u8, shift: bool, no_scale: bool) -> PileUp {
    // read
    let file = File::open(seqfile).ok().expect("read file");
    let mut table = match SeqTable::open(file) {
        Ok(value) => value,
        Err(e) => {
            println!("Error:tabulate: {}", e.to_string()); 
            exit(1); 
        },
    };
    
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
    let minus_shift = if shift {
        let res = (table.params.plus_offset as i16 - (table.params.cut_length as i16 - table.params.minus_offset as i16 - 1i16)) as i32;
        println!("# minus strand shift = {} bp", res);
        res
    } else {
        0i32
    };
    
    let mut iter = bam.records().peekable();
    let mut pileup = PileUp::new(seqinfos, minus_shift, no_scale);
    let scale = compute_scale_factors(&counts);
    
    while pileup.add_data(&mut table, &mut iter, &mut cur_tid, &map, minqual, &scale) {}
    
    pileup
}
