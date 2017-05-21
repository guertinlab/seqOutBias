//!
//!	This module contains C API exported by seqoutbiaslib.
//!
extern crate libc;

use self::libc::{size_t, calloc};
use scale;
use scale::PileUp;
use counts;
use std::slice;
use std::ffi::CStr;
use std::ptr;
use std::mem;
use std::ops::Deref;

#[repr(C)]
pub struct PileUpData(PileUp);

impl Deref for PileUpData {
    type Target = PileUp;
    fn deref(&self) -> &PileUp {
        let PileUpData(ref inner) = *self;
        inner
    }
}

#[repr(C)]
pub struct PileUpPoint {
    plus: f64,
    minus: f64,
}

#[repr(C)]
pub struct Config {
  // add here the flags
  /// Minimum read quality
  min_qual: u8,
  // regions
  // distance range
  // Distance range for included paired reads

  /// Only accept aligned reads that have a mapped pair
  only_paired: bool,
  /// Only accept BAM reads with length equal to 'read-size'
  exact_length: bool,
  /// Use tail edge of reads (3') instead of start edge (5')
  tail_edge: bool,
  /// Shift minus strand counts
  shift_counts: bool,
  /// If false, pileUp represents unscaled counts
  scale_pileup: bool,
}

/// Create a Config structure filled with the same default values as used in the seqOutBias program.
#[no_mangle]
pub extern fn seqoutbias_default_config() -> Config {
  // setup default values for config options
  Config {
    min_qual: 0,
    only_paired: false,
    exact_length: false,
    tail_edge: false,
    shift_counts: false,
    scale_pileup: true,
  }
}

/// Compute genome-wide pile-up for given input seqtbl file and set of BAM files using supplied Config parameters.
#[no_mangle]
pub extern fn seqoutbias_create_pileup(seqtable_filename: *const libc::c_char, bam_filenames: *const *const libc::c_char, n_bams: size_t, config: Config) -> *mut PileUpData {
  let seqtable_filename = unsafe {
    assert!(!seqtable_filename.is_null());
    CStr::from_ptr(seqtable_filename).to_string_lossy().into_owned()
  };

  let bam_filenames = unsafe { 
    assert!(!bam_filenames.is_null());
    slice::from_raw_parts(bam_filenames, n_bams as usize)
  };
  let bams: Vec<String> = bam_filenames.iter()
      .map(|&p| unsafe { CStr::from_ptr(p).to_string_lossy().into_owned() })
      .collect();

  let regions: Option<String> = None;
  let dist_range: Option<(i32, i32)> = None;

  // collect counts
  let counts = counts::tabulate(
    &seqtable_filename, 
    Some(&bams), 
    config.min_qual, 
    regions, 
    dist_range, 
    config.only_paired,
    config.exact_length,
    config.tail_edge
  );

  // compute pileup
  let pileup = scale::scale(
    &seqtable_filename, 
    counts, 
    &bams, 
    config.min_qual,
    config.shift_counts,
    !config.scale_pileup,
    &dist_range,
    config.only_paired,
    config.exact_length,
    config.tail_edge
  );

  Box::into_raw(Box::new(PileUpData(pileup)))
}

/// Free memory used by pile-up
#[no_mangle]
pub extern fn seqoutbias_free_pileup(ptr: *mut PileUpData) {
  if ptr.is_null() { return }
  unsafe { Box::from_raw(ptr); }
}

/// Translate chromosome pile-up into an array of 'PileUpPoint's of length 'out_len'.
/// Caller must free memory using 'free()'.
/// Returns NULL if chromosome name was not found in pile-up.
#[no_mangle]
pub extern fn seqoutbias_chrom_pileup(ptr: *mut PileUpData, chrom: *const libc::c_char, out_len: *mut size_t) -> *mut PileUpPoint {
  let pileup = unsafe {
    assert!(!ptr.is_null());
    &mut *ptr
  };
  let chrom = unsafe {
    assert!(!chrom.is_null());
    CStr::from_ptr(chrom).to_string_lossy()
  };

  // get chrom index
  let chrom_index = match pileup.chrom_index(&chrom) {
    Some(value) => value,
    None => return ptr::null_mut(),
  };

  // get chrom size
  let chrom_size = pileup.chrom_size(chrom_index).unwrap() as usize;

  unsafe {
    assert!(!out_len.is_null());
    *out_len = chrom_size * mem::size_of::<PileUpPoint>();
  }

  // allocate memory with C calloc
  let result: *mut PileUpPoint = unsafe {
    calloc(chrom_size, mem::size_of::<PileUpPoint>()) as *mut PileUpPoint
  };

  let rslice = unsafe {
    slice::from_raw_parts_mut(result, chrom_size)
  };

  // fill memory
  for (pos, value) in pileup.chrom_iter(chrom_index) {
    let idx = *pos as usize;
    rslice[idx].plus = value.0;
    rslice[idx].minus = value.1;
  }
  
  result
}

/// Obtain index (returned via 'out_index') for given chromosome name in pile-up.
/// Returns 1 if chromosome name was found, 0 otherwise.
#[no_mangle]
pub extern fn seqoutbias_query_pileup_chrom_index(ptr: *mut PileUpData, chrom: *const libc::c_char, out_index: *mut size_t) -> i32 {
  let pileup = unsafe {
    assert!(!ptr.is_null());
    &mut *ptr
  };

  let chrom = unsafe {
    assert!(!chrom.is_null());
    CStr::from_ptr(chrom).to_string_lossy()
  };

  match pileup.chrom_index(&chrom) {
    Some(value) => {
      unsafe { 
        assert!(!out_index.is_null());
        *out_index = value; 
      }
      1
    },
    None => 0,
  }
}

/// Query pile-up at given chromosome position (returned via 'out').
/// Chromosome index obtained from 'seqoutbias_query_pileup_chrom_index'.
/// Returns 1 if position has values in pile-up, 0 otherwise ('out' unchanged if return value is 0).
/// NOTE: It's undefined behaviour to pass an invalid chromosome index.
#[no_mangle]
pub extern fn seqoutbias_query_pileup(ptr: *mut PileUpData, chrom_index: size_t, pos: u32, out: *mut PileUpPoint) -> i32 {
  let pileup = unsafe {
    assert!(!ptr.is_null());
    &mut *ptr
  };

  let out = unsafe {
    assert!(!out.is_null());
    &mut *out
  };

  match pileup.get(chrom_index, pos) {
    Some(value) => {
      out.plus = value.0;
      out.minus = value.1;
      1
    },
    None => 0,
  }
}
