extern crate rustc_serialize;
extern crate docopt;
extern crate num;
extern crate flate2;
extern crate byteorder;
extern crate bincode;
extern crate rust_htslib as htslib;

mod tallyrun;
mod tallyread;
mod seqtable;
mod fasta;
mod counts;
mod bigwig;
mod scale;

use docopt::Docopt;
use std::process::exit;
use std::error::Error;

/* Main usage/arguments */

const USAGE: &'static str = "
Cut-site frequencies

Usage:
  enzcut tallymer <fasta-file> <read-size> [--parts=<n>]
  enzcut seqtable <fasta-file> [options]
  enzcut dump <seqtbl-file> [<seqrange>]
  enzcut table <seqtbl-file> [<bam-file>] [--qual=<q>] [--regions=<bedfile>] [--pdist=<min:max>] [--only-paired]
  enzcut scale <seqtbl-file> <bam-file> [options]
  enzcut <fasta-file> [<bam-file>] [options]
  enzcut (-h | --help)
  enzcut --version

Options:
  -h --help             Show this screen.
  --version             Show version.
  --cut-size=<n>        Cut-site size [default: 4].
  --tallymer=<file>     Unmappable positions file produced by tallymer (seq, pos).
  --plus-offset=<p>     Cut-site offset on plus strand [default: 2]. Eg, p=2 AA[A]A.
  --minus-offset=<m>    Cut-site offset on minus strand [default: 2]. Eg, m=2 A[A]AA.
  --read-size=<r>       Read length [default: 36].
  --parts=<n>           Split suffix tree generation into n parts [default: 4].
  --qual=<q>            Minimum read quality [default: 0].
  --regions=<bedfile>   Count only cut-sites inside the regions indicated in the BED file.
  --out=<outfile>       Output seqtable filename [default: output.tbl].
  --bed=<bedfile>       Output scaled BED filename [default: output.bed].
  --bw=<bigwigfile>     Output scaled BigWig filename [default: output.bw].
  --stranded            Output per strand counts when writting scaled values.
  --shift-counts        Shift minus strand counts.
  --no-scale            Skip actual scalling in 'scale' command.
  --pdist=<min:max>     Distance range for included paired reads.
  --only-paired         Only accept aligned reads that have a mapped pair.
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_fasta_file: String,
    arg_read_size: u16,
    arg_seqtbl_file: String,
    arg_seqrange: Option<String>,
    arg_bam_file: Option<String>,
    flag_cut_size: u8,
    flag_tallymer: Option<String>,
    flag_plus_offset: u8,
    flag_minus_offset: u8,
    flag_version: bool,
    flag_read_size: u16,
    flag_parts: u8,
    flag_qual: u8,
    flag_regions: Option<String>,
    flag_out: String,
    flag_stranded: bool,
    flag_bed: String,
    flag_bw: String,
    flag_shift_counts: bool,
    flag_no_scale: bool,
    flag_pdist: Option<String>,
    flag_only_paired: bool,
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
        let counts = counts::tabulate(&args.arg_seqtbl_file, args.arg_bam_file, args.flag_qual, args.flag_regions, dist_range, args.flag_only_paired);
        counts::print_counts(&counts, has_bam);
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
    
    // Process main sequence phases
    
    // phase 1 - tallymer
    let tally_path = if run_tallymer {
        let path = tallyrun::tallymer_createfile(&args.arg_fasta_file, args.flag_read_size, args.flag_parts);
        println!("# tallymer produced/found {:}", path.to_str().unwrap());
        Some(path)
    } else { None };
    
    // phase 2 - seqtable
    let seqtable_file = if run_seqtable {
        let seq_params = seqtable::SeqTableParams {
            cut_length: args.flag_cut_size,
            plus_offset: args.flag_plus_offset,
            minus_offset: args.flag_minus_offset,
            read_length: args.flag_read_size };
        fasta::process_fasta(&args.arg_fasta_file, &tally_path.unwrap(), seq_params, &args.flag_out);
        println!("# seqtable produced {}", &args.flag_out);
        &args.flag_out
    } else {
        &args.arg_seqtbl_file
    };
    
    // phase 3 - tabulate & scale
    if run_scale {
        let bamfile = args.arg_bam_file.clone().unwrap();
        let counts = counts::tabulate(seqtable_file, args.arg_bam_file, args.flag_qual, args.flag_regions, dist_range, args.flag_only_paired);
        let pileup = scale::scale(seqtable_file, counts, bamfile, args.flag_qual, args.flag_shift_counts, args.flag_no_scale);
            
        match pileup.write_bed(&args.flag_bed, args.flag_stranded) {
            Ok(_) => println!("# scale produced {}", &args.flag_bed),
            Err(err) => println!("Error producing BED file: {}", err.description()),
        }
        
        match pileup.write_bw(&args.flag_bw, args.flag_stranded) {
            Ok((f1, f2)) => {
                if f2.is_some() {
                    println!("# scale produced {}", f1);
                    println!("# scale produced {}", f2.unwrap());
                } else {
                    println!("# scale produced {}", &args.flag_bw);    
                }
            },
            Err(err) => println!("Error producing BigWig file: {}", err.description()), 
        }
    }
}
