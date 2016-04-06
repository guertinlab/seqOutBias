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
mod scale;

use docopt::Docopt;
use std::process::exit;

/* Main usage/arguments */

const USAGE: &'static str = "
Cut-site frequencies

Usage:
  enzcut tallymer <fasta-file> <read-size> [--parts=<n>]
  enzcut seqtable <fasta-file> [options]
  enzcut dump <seqtbl-file> [<seqrange>]
  enzcut table <seqtbl-file> [<bam-file>] [--qual=<q>] [--regions=<bedfile>]
  enzcut scale <seqtbl-file> <bam-file> [options]
  enzcut <fasta-file> [<bam-file>] [options]
  enzcut (-h | --help)
  enzcut --version

Options:
  -h --help            Show this screen.
  --version            Show version.
  --cut-size=<n>       Cut-site size [default: 4].
  --tallymer=<file>    Unmappable positions file produced by tallymer (seq, pos).
  --plus-offset=<p>    Cut-site offset on plus strand [default: 2]. Eg, p=2 AA[A]A.
  --minus-offset=<m>   Cut-site offset on minus strand [default: 2]. Eg, m=2 A[A]AA.
  --read-size=<r>      Read length [default: 36].
  --parts=<n>          Split suffix tree generation into n parts [default: 4].
  --qual=<q>           Minimum read quality [default: 0].
  --regions=<bedfile>  Count only cut-sites inside the regions indicated in the BED file.
  --out=<outfile>      Output seqtable filename [default: output.tbl].
  --bed=<bedfile>      Output scaled BED filename [default: output.bed]. 
  --stranded           Output per strand counts when writting scaled values.
  --shift-counts       Shift minus strand counts.
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
    flag_shift_counts: bool,
    cmd_tallymer: bool,
    cmd_seqtable: bool,
    cmd_dump: bool,
    cmd_table: bool,
    cmd_scale: bool,
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
    
    if args.cmd_tallymer {
        let tally_path = tallyrun::tallymer_createfile(&args.arg_fasta_file, args.arg_read_size, args.flag_parts);
        println!("# tallymer produced/found {:}", tally_path.to_str().unwrap());
        
        return;
    }
    
    if args.cmd_seqtable {
        // Find tallymer file (creating it if needed)
        let tally_path = tallyrun::tallymer_createfile(&args.arg_fasta_file, args.flag_read_size, args.flag_parts);
        println!("# tallymer produced/found {:}", tally_path.to_str().unwrap());
        
        // 
        let seq_params = seqtable::SeqTableParams {
        cut_length: args.flag_cut_size,
        plus_offset: args.flag_plus_offset,
        minus_offset: args.flag_minus_offset,
        read_length: args.flag_read_size };

        fasta::process_fasta(&args.arg_fasta_file, &tally_path, seq_params, &args.flag_out);
        println!("# seqtable produced {}", &args.flag_out);
        
        return;
    }
    
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
        let counts = counts::tabulate(&args.arg_seqtbl_file, args.arg_bam_file, args.flag_qual, args.flag_regions);
        counts::print_counts(&counts, has_bam);
        return;
    }
    
    if args.cmd_scale {
        let bamfile = args.arg_bam_file.clone().unwrap();
        let counts = counts::tabulate(&args.arg_seqtbl_file, args.arg_bam_file, args.flag_qual, args.flag_regions);
        let pileup = scale::scale(&args.arg_seqtbl_file, counts, bamfile, args.flag_qual, args.flag_shift_counts);
        
        pileup.write_bed(&args.flag_bed, args.flag_stranded).unwrap();
        println!("# scale produced {}", &args.flag_bed);
        return;
    }
    
    // catch cmd names being interpreted as fasta_file names
    if args.arg_fasta_file.eq("dump") || args.arg_fasta_file.eq("table") || args.arg_fasta_file.eq("tallymer") || args.arg_fasta_file.eq("seqtable") || args.arg_fasta_file.eq("scale") {
        println!("Invalid arguments to {} command.", args.arg_fasta_file);
        println!("{}", USAGE);
        exit(1);
    }
    
    
    //
    // Process FASTA file to create sequence table and the counts table
    //
    
    // Find tallymer file (creating it if needed)
    let tally_path = tallyrun::tallymer_createfile(&args.arg_fasta_file, args.flag_read_size, args.flag_parts);
    println!("# tallymer produced/found {:}", tally_path.to_str().unwrap());
    
    // Parameters
    let seq_params = seqtable::SeqTableParams {
	   cut_length: args.flag_cut_size,
	   plus_offset: args.flag_plus_offset,
	   minus_offset: args.flag_minus_offset,
	   read_length: args.flag_read_size };

    // Generate sequence table
    fasta::process_fasta(&args.arg_fasta_file, &tally_path, seq_params, &args.flag_out);
    println!("# seqtable produced {}", &args.flag_out);
    
    // Generate counts table
    let has_bam = args.arg_bam_file.is_some();
    let counts = counts::tabulate(&args.flag_out, args.arg_bam_file, args.flag_qual, args.flag_regions);
    counts::print_counts(&counts, has_bam);
}
