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
use std::io::BufRead;
use std::io::BufReader;
use std::io::Result;
use std::io::ErrorKind;
use std::io::Error as ioError;
use std::error::Error;
use std::process::exit;
use std::ops::Range;
use std::iter::Peekable;
use seqtable::{SeqTableParams,SeqTable,SequenceInfo};
use std::cmp::Ordering;
use filter::{RecordCheck, PairedChecker, SingleChecker};

struct KeyIter<'a> {
    kmer: Vec<u8>,
    alph: [char; 4],
    mask: Option<&'a Vec<bool>>,
}

impl<'a> KeyIter<'a> {
    fn new(k: u8, mask: Option<&'a Vec<bool>>) -> KeyIter {
        let mut kmer = vec![1; k as usize];
        kmer[(k - 1) as usize] = 0;
        KeyIter { kmer: kmer, alph: ['A','C','G','T'], mask: mask }
    }
    
    fn key(&self) -> String {
        let mut key = String::new();
        
        match self.mask {
            Some(flags) => {
                let kiter = &mut self.kmer.iter();
                for &flag in flags {
                    if flag {
                        let b = kiter.next().unwrap();
                        key.push(self.alph[(b - 1) as usize]);
                    } else {
                        key.push('N');
                    }
                } 
            },
            None => {
                for b in &self.kmer {
                    key.push(self.alph[(b - 1) as usize]);
                }
            }
        }
        key
    }
}

impl<'a> Iterator for KeyIter<'a> {
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
  
  fn contains(&self, chrom_idx: usize, position: i32) -> bool {
      if position < 0 { return false; }
      let position = position as u32;
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
                Some(ref ranges) => ranges.contains(chrom_idx, checker.vir_pos(&record)),
                None => true,
            };
            
            //
            if good && checker.valid(&record) {
                let pair = rdr.vir_get(checker.vir_pos(&record)).unwrap();
                
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
            println!("Error: Failed to read BED file {}: {}", filename, e.description()); 
            exit(1); 
        },
    }
}

fn region_counts<R: ioRead + Seek>(table: &mut SeqTable<R>, bedregions: &str) -> Vec<(u64, u64, u64, u64)> {
    let bediter = bed_regions(bedregions, table);
    
    // allocate counts table
    let mut counts: Vec<(u64, u64, u64, u64)> = Vec::new();
    let nmer_count = table.params.nmer_count();
    
    for _ in 0..nmer_count {
        counts.push((0,0,0,0));
    }
    
    let n_seqs = table.len();
    
    for idx in 0..n_seqs {
        let mut rdr = table.get_sequence_by_idx(idx).ok().expect("read sequence");
        
        for &(start,end) in &bediter.sets[idx] {
            for position in start..end {
                let pair = rdr.get(position).unwrap();
                counts[pair.0 as usize].0 += 1;
                counts[pair.1 as usize].1 += 1;
            }
        }
        
    }
    
    counts
}

fn tabulate_bam<R: ioRead + Seek>(bamfilename: String, seqinfos: &Vec<SequenceInfo>, pair_range: &Option<(i32, i32)>, paired: bool, rlen: usize, minqual: u8, counts: &mut Vec<(u64, u64, u64, u64)>, table: &mut SeqTable<R>, regions: Option<&BedRanges>, exact_length: bool, tail_edge: bool) {
    println!("# tabulate {}", bamfilename);
            
    let bam = match bam::Reader::new(&bamfilename) {
        Ok(value) => value,
        Err(err) => {
            println!("Error: Failed to open BAM '{}': {}", &bamfilename, err.description());
            exit(1);
        },
    };
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
        let checker = match *pair_range {
            Some((min, max)) => PairedChecker {
                tail_edge: tail_edge,
                exact_length: exact_length,
                read_length: rlen,
                min_quality: minqual,
                min_dist: min,
                max_dist: max,
                force_paired: paired,
                max_distance: true,
            },
            None => PairedChecker {
                tail_edge: tail_edge,
                exact_length: exact_length,
                read_length: rlen,
                min_quality: minqual,
                min_dist: 0,
                max_dist: 0,
                force_paired: paired,
                max_distance: false,
            },
        };
        while process_bam_seq(counts, table, &mut iter, &mut cur_tid, &map, &checker, regions) {}
    } else {
        let checker = SingleChecker { tail_edge: tail_edge, exact_length: exact_length, read_length: rlen, min_quality: minqual };
        while process_bam_seq(counts, table, &mut iter, &mut cur_tid, &map, &checker, regions) {}
    }
}

pub fn tabulate(seqfile: &str, bamfile: Option<&Vec<String>>, minqual: u8, regions: Option<String>, pair_range: Option<(i32, i32)>, paired: bool, exact_length: bool, tail_edge: bool) -> Vec<(u64, u64, u64, u64)> {
    // read
    let file = match File::open(seqfile) {
        Ok(value) => value,
        Err(err) => {
            println!("Error: Failed to open sequence table file '{}': {}", seqfile, err.description());
            exit(1);
        },
    };
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
            tabulate_bam(bamfilename.clone(), &seqinfos, &pair_range, paired, rlen, minqual, &mut counts, &mut table, ranges.as_ref(), exact_length, tail_edge);
        }
    }
    
    counts
}

pub fn print_counts(counts: &Vec<(u64, u64, u64, u64)>, with_bam: bool, params: &SeqTableParams) {
    let mut keys = KeyIter::new(params.unmasked_count, params.mask.as_ref());
    
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