use std::fs::File;
use super::read::SeqTable;
use super::read::SequenceInfo;

pub fn dump_seqtable(filename: &str) {
    // read
    let file = File::open(filename).ok().expect("read file");
    let mut table = SeqTable::open(file).ok().expect("read store");
    
    //
    {
        let params = table.params();
        
        println!("# cut-size:     {}", params.cut_length);
        println!("# plus-offset:  {}", params.plus_offset);
        println!("# minus-offset: {}", params.minus_offset);
        println!("# read-size:    {}", params.read_length);
    }
    
    let seqinfos = table.sequences();
    
    for SequenceInfo { name, length } in seqinfos {
        println!(">{} {}", name, length);
        
        let mut rdr = table.get_sequence(&name).ok().expect("read sequence");
        
        for i in 0..8 {
            println!(" {} {:?}", i, rdr.get(i));
        }
        for i in 2000..2008 {
            println!(" {} {:?}", i, rdr.get(i));
        }
    }
}
