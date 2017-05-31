//!
//!	Temporary experimental code.
//!
extern crate rust_htslib as htslib;
use self::htslib::bam;
use self::htslib::bam::Read;
use self::htslib::bam::Reader;
use self::htslib::bam::Records;
use std::collections::BTreeMap;
use std::iter::Peekable;

#[derive(Debug)]
struct SequencePileUp {
    chrom: String,
    counts: BTreeMap<u32, (f64, f64)>
}

impl SequencePileUp {
    // some constructor that takes CountsTable, SeqTable data and BAM files
    // processes the data corresponding to a single sequence and returns a
    // new fully built pileUp

    // op to write pile up into BED file
}


//#[derive(Debug)]
struct BAMSource<'a> {
    bam_set: Vec<Peekable<Records<'a, Reader>>>,
    index: usize,
    id: i32,         // sequence ID inside BAM files ... this may not work ...
    map: Vec<usize>, // map from BAM IDs to SeqTable IDs
}

// The idea of BAMSource is to create a global iterator over Records<Reader>
// that will read all BAM records from a sequence from all input BAM files,
// before moving to the next sequence.
// Want a sort of iterator of iterators

// Basic assumption is that sequences are in order...??


// issues: what happens if not all BAM files have data on all sequences ?
//         do BAM file sequence IDs match across BAM files? probably not!
//         need to create a unified ID table ???

// given a list of input chromosomes:
// - scan BAM headers to map the tid of each one to the chromosome map
// - start with the smallest global id
