extern crate seqoutbiaslib;

use seqoutbiaslib::fasta::reverse_complement;

// Missing data is encoded as zero, so it's reverse complement should also be zero
#[test]
fn kmer_reverse_complement_of_zero_is_zero() {
    assert_eq!( reverse_complement(0, 1), 0 );
}

#[test]
fn kmer_reverse_complement_kmer_size_1() {
    let kmer_size = 1;
    for kmer in 1..5 {
        let rev_comp = 5 - kmer;
        assert_eq!( reverse_complement( kmer, kmer_size), rev_comp );
    }
}

#[test]
fn kmer_reverse_complement_kmer_size_2() {
    let kmer_size = 2;
    for kmer in 1..17 {
        let rev_comp = 17 - kmer;
        assert_eq!( reverse_complement( kmer, kmer_size), rev_comp );
    }
}

#[test]
fn kmer_reverse_complement_kmer_size_4() {
    let kmer_size = 4;
    let kmer = 1u32; // AAAA
    // TTTT
    let rev_comp = 1 + 4u32.pow( 3 ) * 3 + 4u32.pow( 2 ) * 3 + 4u32.pow( 1 ) * 3 + 3;
    assert_eq!( reverse_complement( kmer, kmer_size), rev_comp );
}