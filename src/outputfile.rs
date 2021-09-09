//!
//! Abstraction over output filename composition
//!

use std::path::Path;
use std::ffi::{OsStr, OsString};

#[derive(Debug, Clone)]
pub struct OutFilename {
  stem: OsString,
  extension: OsString,
  suffixes: Vec<OsString>
}

impl OutFilename {

  /// Create new OutFilename instance
  ///
  /// Filename stem will be taken from the `stem_src` (removing the file extension) and the
  /// filename extension will be set to `extension`
  ///
  /// # Arguments
  /// * `stem_src` - file stem source
  /// * `extension` - file extension
  pub fn from_parts( stem_src: &str, extension: &str) -> OutFilename {
    OutFilename {
      stem: Path::new(stem_src).file_stem().unwrap().to_os_string(),
      extension: OsString::from(extension),
      suffixes: Vec::new()
    }
  }

  /// Create new OutFilename instance
  ///
  /// The input filename will be split into a `stem` and `extension`.
  ///
  /// # Arguments
  /// * `filename` - source filename from which to initialize them stem and extension.
  /// * `def_extension` - default extension in case the filename does not have one
  pub fn from_filename( filename: &str, def_extension: &str ) -> OutFilename {
    let def_ext : &OsStr = OsStr::new(def_extension);
    OutFilename {
      stem: Path::new(filename).file_stem().unwrap().to_os_string(),
      extension: Path::new(filename).extension().unwrap_or(def_ext).to_os_string(),
      suffixes: Vec::new()
    }
  }

  /// Create new OutFilename instance
  ///
  /// Helper method that either calls `from_parts()` or `from_filename()` depending on
  /// the presence of the `filename` argument
  ///
  /// # Arguments
  /// * `stem_src` - file stem source
  /// * `filename` - source filename from which to initialize them stem and extension.
  /// * `def_extension` - default extension in case the filename does not have one
  pub fn from( stem: &str, filename: &Option<String>, def_extension: &str) -> OutFilename {
    match filename {
      None => {
        Self::from_parts( stem, def_extension)
      }
      Some(name) => {
        Self::from_filename( name, def_extension )
      }
    }
  }

  pub fn prepend_suffix<'a, T: Into<&'a OsStr>>(&mut self, suffix: T ) {
    self.suffixes.insert( 0, suffix.into().to_os_string() );
  }

  pub fn append_suffix<'a, T: Into<&'a OsStr>>(&mut self, suffix: T ) {
    self.suffixes.push( suffix.into().to_os_string() );
  }

  /// Compute final output filename
  pub fn filename(&self) -> OsString {
    let mut filename = self.stem.to_os_string();
    for suffix in &self.suffixes {
      filename.push( suffix );
    }
    filename.push( "." );
    filename.push( self.extension.to_os_string() );
    return filename;
  }
}

#[cfg(test)]
mod tests {
  use outputfile::OutFilename;
  use std::ffi::{OsStr, OsString};

  #[test]
  fn test_from_parts_with_no_suffixes_added() {
    let outfilename = OutFilename::from_parts( "stem", "ext" );
    assert_eq!( "stem.ext", outfilename.filename() );
  }

  #[test]
  fn test_from_filename_with_no_suffixes_added() {
    let outfilename = OutFilename::from_filename( "stem.ext", "def" );
    assert_eq!( "stem.ext", outfilename.filename() );
  }

  #[test]
  fn test_from_filename_with_default_extension() {
    let outfilename = OutFilename::from_filename( "stem", "def" );
    assert_eq!( "stem.def", outfilename.filename() );
  }

  #[test]
  fn test_from_parts_with_preprended_suffixes() {
    let mut outfilename = OutFilename::from_parts( "stem", "ext" );
    outfilename.prepend_suffix( OsStr::new( "_SUF1" ) );
    outfilename.prepend_suffix( OsStr::new( "_SUF2" ) );
    assert_eq!( "stem_SUF2_SUF1.ext", outfilename.filename() );
  }

  #[test]
  fn test_from_parts_with_appended_suffixes() {
    let mut outfilename = OutFilename::from_parts( "stem", "ext" );
    outfilename.append_suffix( OsStr::new( "_SUF1" ) );
    outfilename.append_suffix( OsStr::new( "_SUF2" ) );
    assert_eq!( "stem_SUF1_SUF2.ext", outfilename.filename() );
  }
}