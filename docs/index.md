Molecular biology enzymes have nucleic acid preferences for their substrates; the preference of an
enzyme is typically dictated by the sequence at or near the active site of the enzyme. This bias may
result in spurious read count patterns when used to interpret high-resolution molecular genomics data.
The seqOutBias program aims to correct this issue by scaling the aligned read counts by the ratio of genome-wide observed read counts to the expected
sequence based counts for each k-mer. The sequence based k-mer counts take into account mappability
at a given read length using Genome Tools' Tallymer program. The seqOutBias program allows for
flexibility in specifying the k-mer, including varying the k-mer size, k-mer information spacing, and
specifying strand-specific offsets for the start of the sequence reads. Due to the large size of some
datasets, seqOutBias reads compressed files (FASTA, mappability information, and BAM files), and
reuses intermediate results as much as possible.

# Requirements

- Platform: OS X or Linux
- Compiler: rust >= 1.5.0 + cargo ( [http://www.rust-lang.org](http://www.rust-lang.org) )
- Genome tools ( [http://genometools.org](http://genometools.org) )
- wigToBigWig ( [UCSC Git Repository](http://genome.ucsc.edu/admin/git.html) or [Precompiled binaries](http://hgdownload.cse.ucsc.edu/admin/exe/))

# Documentation

- [User guide](https://guertinlab.github.io/seqOutBias/seqOutBias_user_guide.pdf)
- [Vignette](https://guertinlab.github.io/seqOutBias/seqOutBias_vignette.pdf)
