//!
//! Code to create a scaled read-count track 
//!
use htslib::bam;
use htslib::bam::Read;
use htslib::bam::Reader;
use htslib::bam::Records;
use std::path::Path;
use std::error::Error;
use std::fs::File;
use std::ffi::OsString;
use std::ffi::OsStr;
use std::io::Read as ioRead;
use std::io::Write as ioWrite;
use std::io::Seek;
use std::io::Error as ioError;
use std::process::exit;
use std::iter::Peekable;
use seqtable::{SeqTable,SequenceInfo};
use std::collections::BTreeMap;
use std::collections::btree_map::Iter;
use bigwig::write_bigwig;
use bigwig::Strand;
use filter::{RecordCheck, PairedChecker, SingleChecker};
use fasta::reverse_complement;

#[derive(Debug)]
pub struct PileUp {
    chroms: Vec<String>,
    chrom_sizes: Vec<u32>,
    counts: Vec<BTreeMap<u32, (f64, f64)>>,
    plus_shift: i32,
    minus_shift: i32,
    no_scale: bool,
}

impl PileUp {
    
    fn new(sinfos: &Vec<SequenceInfo>, plus_shift: i32, minus_shift: i32, no_scale: bool) -> PileUp {
        let mut chroms = Vec::new();
        let mut counts = Vec::new();
        let mut chrom_sizes = Vec::new();
        
        for sinfo in sinfos {
            chroms.push(sinfo.name.clone());
            chrom_sizes.push(sinfo.length);
            counts.push(BTreeMap::new());
        }
        
        PileUp {
            chroms: chroms,
            chrom_sizes: chrom_sizes,
            counts: counts,
            plus_shift: plus_shift,
            minus_shift: minus_shift,
            no_scale: no_scale
        }
    }
    
    fn add_data<R: ioRead+Seek, C: RecordCheck>(&mut self, table: &mut SeqTable<R>, bamrecs: &mut Peekable<Records<Reader>>, tid: &mut i32, map: &Vec<usize>, scale: &Vec<(f64, f64)>, checker: &C) -> bool {
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
        let slen = table.len_by_idx(sidx).ok().expect("read sequence length") as i32;
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
                if checker.valid(&record) {
                    let pos = checker.vir_pos(&record);
                    if pos < slen {
                        let (plus_idx, minus_idx) = rdr.vir_get(pos).unwrap();
                        
                        if record.is_reverse() {
                            if minus_idx == 0 {
                                /* no data */
                            } else {
                                let inc = if self.no_scale { 1f64 } else { scale[minus_idx as usize].1 };
                                let minus_pos = (checker.vir_pos(&record) + rlen as i32 - 1i32 + self.minus_shift) as u32;
                                self.counts[sidx as usize].entry(minus_pos).or_insert((0f64, 0f64)).1 += inc;
                            }
                        } else {
                            if plus_idx == 0 {
                                /* no data */
                            } else {
                                let inc = if self.no_scale { 1f64 } else { scale[plus_idx as usize].0 };
                                let plus_pos = (checker.vir_pos(&record) + self.plus_shift) as u32;
                                self.counts[sidx as usize].entry(plus_pos).or_insert((0f64, 0f64)).0 += inc;
                            }
                        }
                    }
                }
            } 
        }
    }
    
    pub fn write_bed(&self, filename: &str, stranded: bool) -> Result<(),ioError> {
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
    
    fn bw_stranded_filename(filename: &str, strand: Strand) -> OsString {
        let def_ext : &OsStr = OsStr::new("bw");
        let mut basename: OsString = Path::new(filename).file_stem().unwrap_or(OsStr::new(filename)).to_os_string();
        let extension = Path::new(filename).extension().unwrap_or(def_ext).to_os_string();
        
        match strand {
            Strand::Plus => basename.push("_plus."),
            Strand::Minus => basename.push("_minus."),
            Strand::Both => unreachable!(),  
        }
        basename.push(extension);
        
        return basename;
    }
    
    pub fn write_bw(&self, filename: &str, stranded: bool) -> Result<(String
    , Option<String>), ioError> {
        if stranded {
            let output_plus = PileUp::bw_stranded_filename(filename, Strand::Plus);
            let output_minus = PileUp::bw_stranded_filename(filename, Strand::Minus);
            
	        match write_bigwig(&output_plus, &self.chroms, &self.chrom_sizes, &self.counts, Strand::Plus) {
                Ok(_) => {
                    try!(write_bigwig(&output_minus, &self.chroms, &self.chrom_sizes, &self.counts, Strand::Minus));
                    Ok((output_plus.into_string().unwrap(), Some(output_minus.into_string().unwrap())))
                },
                Err(err) => Err(err),
            }
        } else {
            try!(write_bigwig(OsStr::new(filename), &self.chroms, &self.chrom_sizes, &self.counts, Strand::Both));
            Ok((String::from(filename), None))
        }
    }

    pub fn chrom_index(&self, chrom: &str) -> Option<usize> {
        self.chroms.iter().position(|x| x == chrom)
    }

    pub fn chrom_size(&self, index: usize) -> Option<u32> {
        self.chrom_sizes.get(index).map(|&x| x)
    }

    pub fn get(&self, chrom_index: usize, position: u32) -> Option<&(f64, f64)> {
        if chrom_index >= self.counts.len() { return None; }
        self.counts[chrom_index].get(&position)
    }

    pub fn chrom_iter(&self, index: usize) -> Iter<u32, (f64, f64)> {
        self.counts[index].iter()
    }
}

fn scale_factor(exp: u64, etotal: u64, obs: u64, ototal: u64) -> f64 {
    let fexp = if etotal > 0 { exp as f64 / etotal as f64 } else { 0f64 };
    let fobs = if ototal > 0 { obs as f64 / ototal as f64 } else { 0f64 };
    if fobs > 0f64 { fexp / fobs } else { 0f64 }
}

fn compute_scale_factors(counts: &Vec<(u64, u64, u64, u64)>) -> Vec<(f64, f64)> {
    // compute totals
    // skip first row which contains unmappable/unusable position counts
    let totals = counts.iter().skip(1).fold((0u64, 0u64, 0u64, 0u64),
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

pub fn scale(seqfile: &str, counts: Vec<(u64, u64, u64, u64)>, bamfiles: &Vec<String>, minqual: u8, shift: bool, shift_amounts: &Option<(i32, i32)>, no_scale: bool, pair_range: &Option<(i32, i32)>, paired: bool, exact_length: bool, tail_edge: bool) -> PileUp {
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
    
    let seqinfos = table.sequences();
    let rlen = table.params.read_length as usize;
    let minus_shift = if shift {
        let res = (table.params.plus_offset as i16 - (table.params.kmer_length as i16 - table.params.minus_offset as i16 - 1i16)) as i32;
        println!("# minus strand shift = {} bp", res);
        res
    } else {
        if let Some(( _, minus)) = *shift_amounts {
            minus
        } else {
            0i32
        }
    };
    let plus_shift = if let Some((plus, _)) = *shift_amounts {
        plus
    } else {
        0i32
    };

    let mut pileup = PileUp::new(&seqinfos, plus_shift, minus_shift, no_scale);
    
    for bamfile in bamfiles {
        println!("# scale {}", &bamfile);
        let bam = match bam::Reader::from_path(&Path::new(bamfile)) {
            Ok(value) => value,
            Err(err) => {
                println!("Error: Failed to open BAM file '{}': {}", bamfile, err.description());
                exit(1);
            },
        };
        let names = bam.header().target_names();
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
        let scale = compute_scale_factors(&counts);
        
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
            while pileup.add_data(&mut table, &mut iter, &mut cur_tid, &map, &scale, &checker) {}
        } else {
            let checker = SingleChecker { tail_edge: tail_edge, exact_length: exact_length, read_length: rlen, min_quality: minqual };
            while pileup.add_data(&mut table, &mut iter, &mut cur_tid, &map, &scale, &checker) {}
        }
    }
    
    pileup
}
