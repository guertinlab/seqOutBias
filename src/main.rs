#![allow(non_snake_case)]

extern crate seqoutbiaslib;
extern crate rustc_serialize;
extern crate docopt;

use seqoutbiaslib::tallyrun;
use seqoutbiaslib::seqtable;
use seqoutbiaslib::fasta;
use seqoutbiaslib::counts;
use seqoutbiaslib::scale;
use seqoutbiaslib::file_exists;

use docopt::Docopt;
use std::process::exit;
use std::error::Error;
use std::path::Path;
use std::fs::File;
use seqoutbiaslib::seqtable::SeqTable;

/* Main usage/arguments */

const USAGE: &'static str = "
Cut-site frequencies

Usage:
  seqOutBias tallymer <fasta-file> <read-size> [--parts=<n>]
  seqOutBias seqtable <fasta-file> [options]
  seqOutBias dump <seqtbl-file> [<seqrange>]
  seqOutBias table <seqtbl-file> [<bam-file>...] [--qual=<q>] [--regions=<bedfile>] [--pdist=<min:max>] [--only-paired] [--exact-length] [--tail-edge]
  seqOutBias scale <seqtbl-file> <bam-file>... [options]
  seqOutBias <fasta-file> <bam-file>... [options]
  seqOutBias (-h | --help)
  seqOutBias --version

Options:
  -h --help             Show this screen.
  --version             Show version.
  --kmer-size=<n>       Kmer size [default: 4].
  --tallymer=<file>     Unmappable positions file produced by tallymer (seq, pos).
  --plus-offset=<p>     Cut-site offset on plus strand, eg. p=2 AA[A]A [default: 2].
  --minus-offset=<m>    Cut-site offset on minus strand, eg. Eg, m=2 A[A]AA [default: 2].
  --kmer-mask=<str>     String indicating relevant kmer positions and cut-site, eg. NNXXNNCXXXXNNXXNN.
  --read-size=<r>       Read length [default: 36].
  --parts=<n>           Split suffix tree generation into n parts [default: 4].
  --qual=<q>            Minimum read quality [default: 0].
  --regions=<bedfile>   Count only cut-sites inside the regions indicated in the BED file.
  --out=<outfile>       Output seqtable filename (defaults to fasta file basename with .tbl extension).
  --bed=<bedfile>       Output scaled BED filename (defaults to BAM file basename with '_scaled.bed' extension).
  --skip-bed            Skip creating the BED file output.
  --bw=<bigwigfile>     Output scaled BigWig filename (defaults to BAM file basename with .bw extension).
  --skip-bw             Skip creating the BigWig file output.
  --stranded            Output per strand counts when writting scaled values.
  --shift-counts        Shift minus strand counts.
  --no-scale            Skip actual scalling in 'scale' command.
  --pdist=<min:max>     Distance range for included paired reads.
  --only-paired         Only accept aligned reads that have a mapped pair.
  --exact-length        Only accept BAM reads with length equal to 'read-size'.
  --tail-edge           Use tail edge of reads (3') instead of start edge (5').
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_fasta_file: String,
    arg_read_size: u16,
    arg_seqtbl_file: String,
    arg_seqrange: Option<String>,
    arg_bam_file: Option<Vec<String>>,
    flag_kmer_size: u8,
    flag_tallymer: Option<String>,
    flag_plus_offset: u8,
    flag_minus_offset: u8,
    flag_kmer_mask: Option<String>,
    flag_version: bool,
    flag_read_size: u16,
    flag_parts: u8,
    flag_qual: u8,
    flag_regions: Option<String>,
    flag_out: Option<String>,
    flag_stranded: bool,
    flag_bed: Option<String>,
    flag_skip_bed: bool,
    flag_bw: Option<String>,
    flag_skip_bw: bool,
    flag_shift_counts: bool,
    flag_no_scale: bool,
    flag_pdist: Option<String>,
    flag_only_paired: bool,
    flag_exact_length: bool,
    flag_tail_edge: bool,
    cmd_tallymer: bool,
    cmd_seqtable: bool,
    cmd_dump: bool,
    cmd_table: bool,
    cmd_scale: bool,
}

fn parse_range(range: &str) -> (i32, i32) {
    let parts: Vec<&str> = range.split(':').collect();
    
    if parts.len() != 2 {
        println!("Invalid distance range: {}", range);
        exit(1);
    }
    let min_val = parts[0].parse::<i32>();
    let max_val = parts[1].parse::<i32>();
    
    if min_val.is_err() || max_val.is_err() {
        println!("Invalid distance range: {}", range);
        exit(1);
    }
    (min_val.unwrap(),max_val.unwrap())
}

fn stem_filename(stem_src: &str, suffix: &str, out_arg: Option<String>) -> String {
    match out_arg {
        Some(out) => out,
        None => {
            let mut basename = Path::new(stem_src).file_stem().unwrap().to_os_string();
            basename.push(suffix);
            basename.into_string().unwrap()
        }
    }
}

fn validate_mask(mask: &str) {
    if let Err(error) = seqtable::SeqTableParams::validate_mask(mask) {
        println!("{}", error);
        exit(1);
    }
}

fn main() {
    // Parse command line arguments
    let args: Args = Docopt::new(USAGE)
                            .and_then(|d| d.decode())
                            .unwrap_or_else(|e| e.exit());

    if args.flag_version {
        println!("Cut-site frequencies, v{}.{}.{}", 
            env!( "CARGO_PKG_VERSION_MAJOR" ),
            env!( "CARGO_PKG_VERSION_MINOR" ),
            env!( "CARGO_PKG_VERSION_PATCH" ) );
        return;
    }
    
    let dist_range = match args.flag_pdist {
        Some(range) => Some(parse_range(&range)),
        None => None
    };
    
    // Check for data output commands
    
    if args.cmd_dump {
        if let Some(seqrange) = args.arg_seqrange {
            seqtable::dump_seqtable_range(&args.arg_seqtbl_file, &seqrange);   
        } else {
            seqtable::dump_seqtable(&args.arg_seqtbl_file);
        }
        
        return;
    }
    
    if args.cmd_table {
        let has_bam = args.arg_bam_file.is_some();
        let counts = counts::tabulate(&args.arg_seqtbl_file, args.arg_bam_file.as_ref(), args.flag_qual, args.flag_regions, dist_range, args.flag_only_paired, args.flag_exact_length, args.flag_tail_edge);
        let params = seqtable::SeqTableParams::from_file(&args.arg_seqtbl_file);
        counts::print_counts(&counts, has_bam, &params);
        return;
    }
    
    // Check for main sequence commands
    let mut run_tallymer = true;
    let mut run_seqtable = true;
    let mut run_scale = true;
    
    if args.cmd_tallymer {
        // just phase 1
        run_seqtable = false;
        run_scale = false;
    } else if args.cmd_seqtable {
        // phase 1 & 2
        run_scale = false;
    } else if args.cmd_scale {
        // just phase 3
        run_tallymer = false;
        run_seqtable = false;
    } else {
        // all three phases
        
        // catch cmd names being interpreted as fasta_file names
        if args.arg_fasta_file.eq("dump") || args.arg_fasta_file.eq("table") || args.arg_fasta_file.eq("tallymer") || args.arg_fasta_file.eq("seqtable") || args.arg_fasta_file.eq("scale") {
            println!("Invalid arguments to {} command.", args.arg_fasta_file);
            println!("{}", USAGE);
            exit(1);
        }
    }
    
    if args.flag_skip_bed && args.flag_skip_bw {
        println!("--skip-bed and --skip-bw cannot be used together");
        exit(1);
    }
    
    // Process main sequence phases
    
    // phase 1 - tallymer
    let tally_path = if run_tallymer {
        if !file_exists(&args.arg_fasta_file) {
            println!("Error: FASTA file {} does not exist!", args.arg_fasta_file);
            exit(1);
        }
        
        let path = tallyrun::tallymer_createfile(&args.arg_fasta_file, args.flag_read_size, args.flag_parts);
        println!("# tallymer produced/found {:}", path.to_str().unwrap());
        Some(path)
    } else { None };
    
    // phase 2 - seqtable
    let seqtable_file = if run_seqtable {

        if let Some(ref mask) = args.flag_kmer_mask {
            validate_mask(mask);
        }

        let seq_params = seqtable::SeqTableParams::new(
            args.flag_kmer_size,
            args.flag_plus_offset,
            args.flag_minus_offset,
            args.flag_read_size,
            &args.flag_kmer_mask);
        
        println!("# kmer-size: {}", seq_params.kmer_length);
        println!("# plus-offset: {}", seq_params.plus_offset);
        println!("# minus-offset: {}", seq_params.minus_offset);
        
        let suffix = format!("_{}.{}.{}.{}.tbl", seq_params.read_length, seq_params.kmer_length, seq_params.plus_offset, seq_params.minus_offset);
        let outfile = stem_filename(&args.arg_fasta_file, &suffix, args.flag_out);
        
        if file_exists(&outfile) {
            let file = File::open(&outfile).ok().expect("read file");
            let table = match SeqTable::open(file) {
                Ok(value) => value,
                Err(e) => {
                    println!("Error: seqtable: {}", e.to_string()); 
                    exit(1);
                },
            };
            if table.equivalent(&outfile, &seq_params) {
                println!("# seqtable reusing existing {}", &outfile);
                outfile
            } else {
                println!("Error: seqtable: output file {} already exists but does not match requested parameters!", outfile); 
                exit(1);
            }
        } else {
            fasta::process_fasta(&args.arg_fasta_file, &tally_path.unwrap(), &seq_params, &outfile);
            println!("# seqtable produced {}", &outfile);
            outfile
        }
    } else {
        args.arg_seqtbl_file
    };
    
    // phase 3 - tabulate & scale
    if run_scale {
        for filename in args.arg_bam_file.as_ref().unwrap() {
            if !file_exists(filename) {
                println!("Error: BAM file {} does not exist!", filename);
                exit(1);
            }
        }
        
        let bamfile = args.arg_bam_file.as_ref().unwrap()[0].clone(); // use the first name for reference
        let counts = counts::tabulate(&seqtable_file, args.arg_bam_file.as_ref(), args.flag_qual, args.flag_regions, dist_range, args.flag_only_paired, args.flag_exact_length, args.flag_tail_edge);
        let pileup = scale::scale(&seqtable_file, counts, args.arg_bam_file.as_ref().unwrap(), args.flag_qual, args.flag_shift_counts, args.flag_no_scale, &dist_range, args.flag_only_paired, args.flag_exact_length, args.flag_tail_edge);
        
        if !args.flag_skip_bed {
            let outfile_bed = if args.flag_no_scale {
                    stem_filename(&bamfile, "_not_scaled.bed", args.flag_bed)
                } else {
                    stem_filename(&bamfile, "_scaled.bed", args.flag_bed)
                };
            
            if file_exists(&outfile_bed) {
                println!("Error: output BED file {} already exists!", outfile_bed);
                exit(1);
            }
            
            match pileup.write_bed(&outfile_bed, args.flag_stranded) {
                Ok(_) => println!("# scale produced {}", &outfile_bed),
                Err(err) => println!("Error producing BED file: {}", err.description()),
            }
        } else {
            println!("# scale skipping BED output");
        }
        
        if !args.flag_skip_bw {
            let outfile_bw = stem_filename(&bamfile, ".bw", args.flag_bw);
            
            match pileup.write_bw(&outfile_bw, args.flag_stranded) {
                Ok((f1, f2)) => {
                    if f2.is_some() {
                        println!("# scale produced {}", f1);
                        println!("# scale produced {}", f2.unwrap());
                    } else {
                        println!("# scale produced {}", f1);   
                    }
                },
                Err(err) => println!("Error producing BigWig file: {}", err.description()), 
            }
        } else {
            println!("# scale skipping BigWig output");
        }
    }
}
