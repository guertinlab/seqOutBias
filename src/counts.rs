//!
//!	This module contains the code to sum up n-mer aligned read counts based on a
//! SeqTable instance. 
//!
use htslib::bam;
use htslib::bam::Read;
use htslib::bam::Reader;
use htslib::bam::Records;
use std::fs::File;
use std::io::Read as ioRead;
use std::io::Seek;
use std::process::exit;
use std::iter::Peekable;
use seqtable::SeqTable;

struct KeyIter {
    kmer: Vec<u8>,
    alph: [char; 4],
}

impl KeyIter {
    fn new(k: u8) -> KeyIter {
        let mut kmer = vec![1; k as usize];
        kmer[(k - 1) as usize] = 0;
        KeyIter { kmer: kmer, alph: ['A','C','G','T'] }
    }
    
    fn key(&self) -> String {
        let mut key = String::new();
        
        for b in &self.kmer {
            key.push(self.alph[(b - 1) as usize]);
        }
        key
    }
}

impl Iterator for KeyIter {
    type Item = String;
    
    fn next(&mut self) -> Option<String> {
        let mut idx = self.kmer.len() - 1;
        
        loop {
            self.kmer[idx] += 1;
            if self.kmer[idx] > 4 {
                self.kmer[idx] = 1;
                
                if idx > 0 {
                  idx -= 1;  
                } else {
                  return None;  
                }
            } else {
                break;
            }
        }
        
        Some(self.key())
    }
}

fn process_bam_seq<R: ioRead+Seek>(counts: &mut Vec<(u64, u64, u64, u64)>, table: &mut SeqTable<R>, bamrecs: &mut Peekable<Records<Reader>>, tid: &mut i32, map: &Vec<usize>, rlen: usize, minqual: u8) -> bool {
    // skip unmapped sequences (tid = -1)
    if *tid < 0 {
        match bamrecs.next() {
            Some(Ok(ref rec)) => { *tid = rec.tid(); return true; },
            Some(Err(_)) => return false,
            None => return false,
        }
    }
    // TODO: store counts table in seqtable file ...
    
    let mut rdr = table.get_sequence_by_idx(map[*tid as usize]).ok().expect("read sequence");
    
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
                let pair = rdr.get(record.pos() as u32).unwrap();
                
                if record.is_reverse() {
                    counts[pair.1 as usize].3 += 1;
                } else {
                    counts[pair.0 as usize].2 += 1;
                } 
            }
        } 
    }
}

pub fn tabulate(seqfile: &str, bamfile: Option<String>, minqual: u8) -> Vec<(u64, u64, u64, u64)> {
    // read
    let file = File::open(seqfile).ok().expect("read file");
    let mut table = match SeqTable::open(file) {
        Ok(value) => value,
        Err(e) => {
            println!("Error:tabulate: {}", e.to_string()); 
            exit(1); 
        },
    };
    
    // allocate counts table
    let mut counts: Vec<(u64, u64, u64, u64)> = Vec::new();
    let nmer_count = 4u64.pow(table.params.cut_length as u32) + 1;
    
    for _ in 0..nmer_count {
        counts.push((0,0,0,0));
    }
    
    let rlen = table.params.read_length as usize;
    let seqinfos = table.sequences();
    
    for idx in 0..seqinfos.len() {
        let length = seqinfos[idx].length;
        let mut rdr = table.get_sequence_by_idx(idx).ok().expect("read sequence");
        
        for i in 0..length {
            let pair = rdr.get(i).unwrap();
            
            counts[pair.0 as usize].0 += 1;
            counts[pair.1 as usize].1 += 1;
        }
    }
    
    // if we received a BAM file, parse it
    if let Some(bamfilename) = bamfile {
        let bam = bam::Reader::new(&bamfilename).ok().expect("Error opening bam.");
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
        
        while process_bam_seq(&mut counts, &mut table, &mut iter, &mut cur_tid, &map, rlen, minqual) {}
    }
    
    counts
}

pub fn print_counts(counts: &Vec<(u64, u64, u64, u64)>, with_bam: bool) {
    let k = (((counts.len() - 1) as f32).log2() / 2.0) as u8;
    let mut keys = KeyIter::new(k);
    
    // TODO: tmp: output table to file
    if with_bam {
        for i in 1..counts.len() {
            let (plus, minus, bam_plus, bam_minus) = counts[i];
            let key = keys.next().unwrap();
            println!("{}\t{}\t{}\t{}\t{}\t{}", i, key, plus, minus, bam_plus, bam_minus);
        }
    } else {
        for i in 1..counts.len() {
            let (plus, minus, _, _) = counts[i];
            let key = keys.next().unwrap();
            println!("{}\t{}\t{}\t{}", i, key, plus, minus);
        }
    }
}