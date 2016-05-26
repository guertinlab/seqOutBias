//!
//!	This module contains the code to sum up n-mer aligned read counts based on a
//! SeqTable instance. 
//!
use htslib::bam;
use htslib::bam::Read;
use htslib::bam::Reader;
use htslib::bam::Records;
use htslib::bam::record::Record;
use std::fs::File;
use std::io::Read as ioRead;
use std::io::Seek;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Result;
use std::io::ErrorKind;
use std::io::Error as ioError;
use std::error::Error;
use std::ops::Range;
use std::process::exit;
use std::iter::Peekable;
use seqtable::SeqTable;
use std::cmp::Ordering;

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


// Data for iterator over BED regions
struct BedRanges {
  sets: Vec<Vec<(u32, u32)>>,
  chrom_idx: usize,
  row_idx: usize,
}

impl BedRanges {
  
  pub fn parse(filename: &str, chroms: Vec<String>) -> Result<BedRanges> {
    let file = try!(File::open(filename));
    let reader = BufReader::new(file);
    let mut n_lines = 0;
    let mut n_no_chrom = 0;
    let mut sets = Vec::new();
    
    for _ in 0..chroms.len() {
      sets.push(Vec::new());
    }
    
    for res in reader.lines() {
      let line = try!(res);
      let parts: Vec<&str> = line.split('\t').collect();
      
      //
      if parts.len() < 3 {
        return Err(ioError::new(ErrorKind::Other, format!("not enough fields in line '{}'", line)));
      }
      n_lines += 1;
      
      // locate chromosome in set
      let idx = match chroms.iter().position(|chr| chr == parts[0]) {
        Some(idx) => idx,
        None => {
          n_no_chrom += 1;
          continue;
        }
      };
      
      // add region to set
      let start: u32 = match parts[1].parse() {
        Ok(value) => value,
        Err(_) => return Err(ioError::new(ErrorKind::Other, format!("invalid field value in line '{}'", line))),
      };
      let end: u32 = match parts[2].parse() {
        Ok(value) => value,
        Err(_) => return Err(ioError::new(ErrorKind::Other, format!("invalid field value in line '{}'", line))),
      };
      sets[idx].push((start, end));
    }
    
    // sort sets by starting position
    for i in 0..sets.len() {
      sets[i].sort_by(|a, b| a.0.cmp(&b.0)); // TODO: consider converting this to sort_by_key - requires rust 1.7
    }
    
    // check that no stored range overlaps
    for i in 0..sets.len() {
      for j in 1..sets[i].len() {
        if sets[i][j].0 < sets[i][j - 1].1 {
          return Err(ioError::new(ErrorKind::Other, format!("overlapping ranges in {}: {}-{} & {}-{}", chroms[i], sets[i][j - 1].0, sets[i][j - 1].1, sets[i][j].0, sets[i][j].1)));
        }
      }
    }
    
    if n_no_chrom > 0 {
      println!("WARN:{}: {} rows ({:.2} %) from unknown chromosomes.", filename, n_no_chrom, n_no_chrom as f32 / n_lines as f32 * 100.0);
    }
    
    Ok(BedRanges{ sets: sets, chrom_idx: 0, row_idx: 0})
  }
  
  fn contains(&self, chrom_idx: usize, position: u32) -> bool {
      self.sets[chrom_idx].binary_search_by(|probe| {
          if probe.0 <= position && probe.1 > position {
              Ordering::Equal
          } else if probe.0 > position {
              Ordering::Greater
          } else {
              Ordering::Less
          }
      }).is_ok()
  }
}

// Iterator over BED regions returns chromosome index (as per supplied map) and the region coordinates in the form of a iteratable Range
impl Iterator for BedRanges {
  type Item = (usize, Range<u32>);
  
  fn next(&mut self) -> Option<(usize, Range<u32>)> {
    if self.row_idx >= self.sets[self.chrom_idx].len() {
      self.row_idx = 0;
      self.chrom_idx += 1;
    }
    if self.chrom_idx >= self.sets.len() {
      return None;
    }
    
    let result = Some((self.chrom_idx, Range { start: self.sets[self.chrom_idx][self.row_idx].0, end: self.sets[self.chrom_idx][self.row_idx].1}));
    self.row_idx += 1;
    result
  }
}

trait RecordCheck {
    fn valid(&self, rec: &Record) -> bool;
}

struct SingleChecker {
    read_length: usize,
    min_quality: u8,
}

impl RecordCheck for SingleChecker {
    fn valid(&self, record: &Record) -> bool {
        !record.is_unmapped() && record.seq().len() == self.read_length && record.mapq() >= self.min_quality
    }
}

struct PairedChecker {
    read_length: usize,
    min_quality: u8,
    min_dist: i32,
    max_dist: i32,
    force_paired: bool,
}

impl RecordCheck for PairedChecker {
    fn valid(&self, record: &Record) -> bool {
        // check single read conditions
        if record.is_unmapped() || record.seq().len() == self.read_length || record.mapq() >= self.min_quality {
            return false;
        }
        // mandatory paired condition
        if ( !record.is_paired() || record.is_mate_unmapped() || record.tid() != record.mtid() ) && self.force_paired {
            return false;
        }
        // check pair distance 
        if record.is_paired() {
            let dist = (record.pos() - record.mpos()).abs() + self.read_length as i32;
            if dist < self.min_dist || dist > self.max_dist {
                return false;
            }
        }
        true
    }
}

fn process_bam_seq<R: ioRead+Seek, C: RecordCheck>(counts: &mut Vec<(u64, u64, u64, u64)>, table: &mut SeqTable<R>, bamrecs: &mut Peekable<Records<Reader>>, tid: &mut i32, map: &Vec<usize>, checker: &C, regions: Option<&BedRanges>) -> bool {
    // skip unmapped sequences (tid = -1)
    if *tid < 0 {
        match bamrecs.next() {
            Some(Ok(ref rec)) => { *tid = rec.tid(); return true; },
            Some(Err(_)) => return false,
            None => return false,
        }
    }
    
    let chrom_idx = map[*tid as usize];
    let mut rdr = table.get_sequence_by_idx(chrom_idx).ok().expect("read sequence");
    
    loop {
        // check if we changed sequence
        match bamrecs.peek() {
            Some(&Ok(ref rec)) => if rec.tid() != *tid { *tid = rec.tid(); return true; },
            Some(&Err(_)) => return false,
            None => return false,
        }
        
        // if not count position
        if let Some(Ok(record)) = bamrecs.next() {
            // point in regions
            let good = match regions.as_ref() {
                Some(ref ranges) => ranges.contains(chrom_idx, record.pos() as u32),
                None => true,
            };
            
            //
            if good && checker.valid(&record) {
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

fn bed_regions<R: ioRead + Seek>(filename: &str, table: &mut SeqTable<R>) -> BedRanges {
    let chroms: Vec<String> = table.sequences().iter().map(|sinfo| sinfo.name.clone()).collect();
    
    match BedRanges::parse(filename, chroms) {
        Ok(iter) => iter,
        Err(e) => {
            println!("Error reading BED file {}: {}", filename, e.description()); 
            exit(1); 
        },
    }
}

fn region_counts<R: ioRead + Seek>(table: &mut SeqTable<R>, bedregions: &str) -> Vec<(u64, u64, u64, u64)> {
    let bediter = bed_regions(bedregions, table);
    
    // allocate counts table
    let mut counts: Vec<(u64, u64, u64, u64)> = Vec::new();
    let nmer_count = 4u64.pow(table.params.cut_length as u32) + 1;
    
    for _ in 0..nmer_count {
        counts.push((0,0,0,0));
    }
    
    for (idx, range) in bediter {
        let mut rdr = table.get_sequence_by_idx(idx).ok().expect("read sequence");
        
        for position in range {
            let pair = rdr.get(position).unwrap();
            counts[pair.0 as usize].0 += 1;
            counts[pair.1 as usize].1 += 1;
        } 
    }
    
    counts
}

pub fn tabulate(seqfile: &str, bamfile: Option<&Vec<String>>, minqual: u8, regions: Option<String>, pair_range: Option<(i32, i32)>, paired: bool) -> Vec<(u64, u64, u64, u64)> {
    // read
    let file = File::open(seqfile).ok().expect("read file");
    let mut table = match SeqTable::open(file) {
        Ok(value) => value,
        Err(e) => {
            println!("Error:tabulate: {}", e.to_string()); 
            exit(1); 
        },
    };
    
    // get counts table from file
    let mut counts = match regions.as_ref() {
        Some(regfile) => region_counts(&mut table, regfile),
        None => table.counts().unwrap(),  
    };
    
    //
    let rlen = table.params.read_length as usize;
    let seqinfos = table.sequences();
        
    // if we received a BAM file, parse it
    if let Some(bamfilenames) = bamfile {
        // bed regions
        let ranges = if let Some(regfile) = regions {
            Some(bed_regions(&regfile, &mut table))
        } else {
            None
        };
        
        for bamfilename in bamfilenames {
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
            
            if pair_range.is_some() || paired {
                let checker = PairedChecker {
                    read_length: rlen,
                    min_quality: minqual,
                    min_dist: 0,
                    max_dist: 0,
                    force_paired: paired,
                };
                while process_bam_seq(&mut counts, &mut table, &mut iter, &mut cur_tid, &map, &checker, ranges.as_ref()) {}
            } else {
                let checker = SingleChecker { read_length: rlen, min_quality: minqual };
                while process_bam_seq(&mut counts, &mut table, &mut iter, &mut cur_tid, &map, &checker, ranges.as_ref()) {}
            }
        }
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