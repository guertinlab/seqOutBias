extern crate rustc_serialize;
extern crate docopt;
extern crate num;
extern crate flate2;
extern crate byteorder;
extern crate bincode;

mod tallyrun;
mod tallyread;
mod seqtable;
mod fasta;

use docopt::Docopt;

/* Main usage/arguments */

const USAGE: &'static str = "
Cut-site frequencies

Usage:
  enzcut tallymer <fasta-file> <read-size> [--parts=<n>]
  enzcut dump <seqtbl-file>
  enzcut <fasta-file> [options]
  enzcut (-h | --help)
  enzcut --version

Options:
  -h --help           Show this screen.
  --version           Show version.
  --cut-size=<n>      Cut-site size [default: 4].
  --tallymer=<file>   Unmappable positions file produced by tallymer (seq, pos).
  --plus-offset=<p>   Cut-site offset on plus strand [default: 2]. Eg, p=2 AA[A]A.
  --minus-offset=<m>  Cut-site offset on minus strand [default: 2]. Eg, m=2 A[A]AA.
  --read-size=<r>     Read length [default: 36].
  --parts=<n>         Split suffix tree generation into n parts [default: 4].
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_fasta_file: String,
    arg_read_size: u16,
    arg_seqtbl_file: String,
    flag_cut_size: u8,
    flag_tallymer: Option<String>,
    flag_plus_offset: u8,
    flag_minus_offset: u8,
    flag_version: bool,
    flag_read_size: u16,
    flag_parts: u8,
    cmd_tallymer: bool,
    cmd_dump: bool,
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
        println!("produced/found {:}", tally_path.to_str().unwrap());
        
        return;
    }
    
    if args.cmd_dump {
        seqtable::dump_seqtable(&args.arg_seqtbl_file);
        
        return;
    }
    
    //
    // Process FASTA file to create sequence table
    //
    
    // Find tallymer file (creating it if needed)
    let tally_path = tallyrun::tallymer_createfile(&args.arg_fasta_file, args.flag_read_size, args.flag_parts);
    println!("produced/found {:}", tally_path.to_str().unwrap());
    
    // 
    let seq_params = seqtable::SeqTableParams {
	   cut_length: args.flag_cut_size,
	   plus_offset: args.flag_plus_offset,
	   minus_offset: args.flag_minus_offset,
	   read_length: args.flag_read_size };

    fasta::process_fasta(&args.arg_fasta_file, &tally_path, seq_params, "output.tbl");
}
