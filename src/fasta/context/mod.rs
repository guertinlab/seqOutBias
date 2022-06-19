//!
//! Code to track current kmer index as new DNA bases are observed
//!

mod simple;
mod masked;
mod masked_stranded;

pub use self::simple::EnzContextSimple;
pub use self::masked::EnzContextMasked;
pub use self::masked_stranded::EnzContextMaskedStrandSpecific;

/// Observed kmer on the plus and minus strands
/// Unless the strand-specific flag is active, both values are equal
/// and match the plus strand readout.
pub struct KmerIndex {
    pub plus: Option<u32>,
    pub minus: Option<u32>
}

/// Holds information on the current enzyme cut region
pub trait EnzContext {
    fn sequence_change(&mut self);
    fn add_base(&mut self, base: u8) -> KmerIndex;
}
