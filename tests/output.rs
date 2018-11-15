extern crate assert_cmd;
extern crate tempdir;
extern crate fs_extra;

use fs_extra::dir::copy;
use fs_extra::dir::CopyOptions;
use std::path::PathBuf;
use std::path::Path;
use std::fs::{self, DirEntry};
use assert_cmd::prelude::*;
use std::process::Command;
use tempdir::TempDir;

fn get_resource_folder() -> PathBuf {
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests");
    d.push("resources");
    return d;
}

fn setup_test_folder() -> TempDir {
    let res = get_resource_folder();
    let parent = res.parent().unwrap();
    let tmp_dir = TempDir::new_in(parent, "test").unwrap();
    
    {
        let tmp_path = tmp_dir.path();

        let options = CopyOptions::new(); //Initialize default values for CopyOptions
        copy(&res, tmp_path, &options).unwrap();
    }

    return tmp_dir;
}

#[test]
fn temporary_bw_files_are_removed() {
    let aux = setup_test_folder();
    let mut src_path = aux.path().to_path_buf();
    src_path.push("resources");

    println!("Test files in {:?}", src_path);

    let mut cmd = Command::main_binary().unwrap();
    cmd.current_dir(&src_path)
       .arg("ref2.fa")
       .arg("reads.bam")
       .arg("--read-size=5")
       .arg("--skip-bed")
       .assert().success();
    
    
    // assert some files exist
    let mut bwOut = src_path.clone();
    bwOut.push("reads.bw");
    assert!(bwOut.exists());
    // reads.bw

    // assert some files are removed
    // chromInfo*.tmp
    // wiggle*.tmp
    for entry in fs::read_dir(&src_path).unwrap() {
        let entry = entry.unwrap();
        let filename = entry.file_name().into_string().unwrap();

        if ( filename.starts_with("chromInfo") || filename.starts_with("wiggle") ) && filename.ends_with(".tmp") {
            panic!("Temporary file has not been removed: {}", filename);
        } 
    }
}
