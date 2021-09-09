extern crate assert_cmd;
extern crate predicates;
extern crate tempdir;
extern crate fs_extra;

use fs_extra::dir::copy;
use fs_extra::dir::CopyOptions;
use fs_extra::file::move_file;
use std::path::PathBuf;
use std::fs;
use assert_cmd::prelude::*;
use std::process::Command;
use tempdir::TempDir;

fn get_test_folder() -> PathBuf {
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests");
    return d;
}

fn get_resource_folder() -> PathBuf {
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests");
    d.push("resources");
    return d;
}

fn setup_test_folder() -> TempDir {
    let mut res = get_resource_folder();
    res.push("base");
    let parent = get_test_folder();
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
    src_path.push("base");

    println!("Test files in {:?}", src_path);

    let mut cmd = Command::main_binary().unwrap();
    cmd.current_dir(&src_path)
       .arg("ref2.fa")
       .arg("reads.bam")
       .arg("--read-size=5")
       .arg("--skip-bed")
       .assert().success();
    
    
    // assert some files exist
    // reads.bw
    let mut bw_out = src_path.clone();
    bw_out.push("reads.bigWig");
    assert!(bw_out.exists());

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

#[test]
fn supplied_tallymer_path_is_respected() {
    let aux = setup_test_folder();
    let mut src_path = aux.path().to_path_buf();
    src_path.push("base");

    // get tallymer path
    let mut tallymer_path = get_resource_folder();
    tallymer_path.push("tallymer");
    tallymer_path.push("ref2.tal_5.gtTxt.gz");

    println!("Test files in {:?}", src_path);

    let mut cmd = Command::main_binary().unwrap();
    cmd.current_dir(&src_path)
       .arg("ref2.fa")
       .arg("reads.bam")
       .arg("--read-size=5")
       .arg("--skip-bed")
       .arg(format!("--tallymer={}", tallymer_path.to_str().unwrap()))
       .assert().success();
    
    // check that tallymer file was not created
    let mut tally_bad = src_path.clone();
    tally_bad.push("ref2.tal_5.gtTxt.gz");
    assert!(!tally_bad.exists());

    // assert some files exist
    // reads.bw
    let mut bw_out = src_path.clone();
    bw_out.push("reads.bigWig");
    assert!(bw_out.exists());
}

fn file_exists_in_folder(filename: &str, path: &PathBuf) -> bool {
    let mut file_path = path.clone();
    file_path.push(filename);
    return file_path.exists();
}

#[test]
fn user_specified_genome_tools_working_directory() {
    // setup two folders
    // 1. folder containing refseq dataset and where genome tools will operate
    let mut res = get_resource_folder();
    res.push("base");
    let parent = get_test_folder();
    let gt_folder = TempDir::new_in(&parent, "test_gt").unwrap();
    
    {
        let tmp_path = gt_folder.path();

        let options = CopyOptions::new(); //Initialize default values for CopyOptions
        copy(&res, tmp_path, &options).unwrap();
    }
    let mut gt_path = gt_folder.path().to_path_buf();
    gt_path.push("base");

    // 2. folder where main program will operate
    let main_folder = TempDir::new_in(&parent, "test_main").unwrap();
    let main_path = main_folder.path().to_path_buf();
    let move_options = fs_extra::file::CopyOptions::new();
    let mut bam_src = gt_path.clone();
    bam_src.push("reads.bam");
    let mut bam_tgt = main_path.clone();
    bam_tgt.push("reads.bam");
    move_file(bam_src, bam_tgt, &move_options).unwrap();

    // execute command
    let mut ref_path = gt_path.clone();
    ref_path.push("ref2.fa");

    let mut cmd = Command::main_binary().unwrap();
    cmd.current_dir(&main_path)
       .arg(ref_path)
       .arg("reads.bam")
       .arg("--read-size=5")
       .arg("--skip-bed")
       .arg(format!("--gt-workdir={}", gt_path.to_str().unwrap()))
       .assert().success();
    
    // check that output files exist in the main directory
    assert!(file_exists_in_folder("reads.bigWig", &main_path));
    assert!(file_exists_in_folder("ref2.tal_5.gtTxt.gz", &main_path));

    // check that genome tools temporary files do not exist in the main directory
    assert!(!file_exists_in_folder("ref2.fa.sa", &main_path));
    assert!(!file_exists_in_folder("ref2.sft.suf", &main_path));
}
