extern crate tempfile;

use std::collections::BTreeMap;
use std::io::Error;
use std::io::ErrorKind;
use std::io::Write;
use std::ffi::OsStr;
use std::fs::File;
use std::fs::remove_file;
use std::process::Command;
//use self::tempfile::Builder;

use ::randfile::random_file;

pub enum Strand {
    Plus,
    Minus,
    Both
}

fn write_chrom_info(filename: &str, chroms: &Vec<String>, chrom_sizes: &Vec<u32>) -> Result<(), Error> {
    let mut f = try!(File::create(filename));
    
    for index in 0..chroms.len() {
        try!(write!(f, "{}\t{}\n", chroms[index], chrom_sizes[index]));
    }
    
    Ok(())
}

fn write_wiggle(filename: &str, chroms: &Vec<String>, counts: &Vec<BTreeMap<u32, (f64, f64)>>, strand: Strand) -> Result<(), Error> {
    // open new file
    let mut f = try!(File::create(filename));
    
    // write data
    for i in 0..chroms.len() {
        let chrom : &str = &chroms[i];
        
        try!(write!(f, "variableStep chrom={}\n", chrom));
        
        for (pos, value) in counts[i].iter() {
            let pos = pos + 1; // shift position by one since wiggleFiles are use 1 based positions
            
            match strand {
                Strand::Plus => if value.0 > 0f64 { try!(write!(f, "{}\t{}\n", pos, value.0)) },
                Strand::Minus => if value.1 > 0f64 { try!(write!(f, "{}\t{}\n", pos, -value.1)) },
                Strand::Both => try!(write!(f, "{}\t{}\n", pos, value.0 + value.1)),
            }
        }
    }
    Ok(())
}

pub fn write_bigwig(filename: &OsStr, chroms: &Vec<String>, chrom_sizes: &Vec<u32>, counts: &Vec<BTreeMap<u32, (f64, f64)>>, strand: Strand) -> Result<(),Error> {
    // create temporary chromInfo file
    /*let chrominfo_tmpfile = Builder::new()
                        .prefix("chromInfo")
                        .suffix(".tmp")
                        .rand_bytes(5)
                        .tempfile()
                        .unwrap();*/
    let chrominfo_tmpfile = random_file("chromInfo", ".tmp", 5);
    
    // create temporary wiggle file
    /*let wiggle_tmpfile = Builder::new()
                        .prefix("wiggle")
                        .suffix(".tmp")
                        .rand_bytes(5)
                        .tempfile()
                        .unwrap();*/
    let wiggle_tmpfile = random_file("wiggle", ".tmp", 5);

    { // block constrains lifetime of borrow on temporary file names
        // write to temporary files
        let chrominfo_name = chrominfo_tmpfile.as_path()//.path()
                            .file_name().unwrap()
                            .to_str().unwrap();
        try!(write_chrom_info(chrominfo_name, chroms, chrom_sizes));

        let wiggle_name = wiggle_tmpfile.as_path()//.path()
                        .file_name().unwrap()
                        .to_str().unwrap();
        try!(write_wiggle(wiggle_name, chroms, counts, strand));
        
        // run wigToBigWig
        let status = Command::new("wigToBigWig")
            .arg("-keepAllChromosomes")
            .arg("-clip")
            .arg(wiggle_name)
            .arg(chrominfo_name)
            .arg(filename)
            .status().unwrap_or_else(|e| {
                panic!("failed to execute wigToBigWig: {}", e)
        });
        
        //
        if !status.success() {
            return Err(Error::new(ErrorKind::Other, format!("An error occurred running: wigToBigWig -keepAllChromosomes {:?} {:?} {:?}", filename, wiggle_name, chrominfo_name)));
        }
    }

    // force release of temporary files
    //wiggle_tmpfile.close().unwrap();
    //chrominfo_tmpfile.close().unwrap();
    remove_file(wiggle_tmpfile.as_path()).expect(&format!("Failed to delete {:?}", wiggle_tmpfile));
    remove_file(chrominfo_tmpfile.as_path()).expect(&format!("Failed to delete {:?}", chrominfo_tmpfile));
    Ok(())
}
