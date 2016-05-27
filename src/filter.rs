use htslib::bam::record::Record;

pub trait RecordCheck {
    fn valid(&self, rec: &Record) -> bool;
}

pub struct SingleChecker {
    pub read_length: usize,
    pub min_quality: u8,
}

impl RecordCheck for SingleChecker {
    fn valid(&self, record: &Record) -> bool {
        !record.is_unmapped() && record.seq().len() == self.read_length && record.mapq() >= self.min_quality
    }
}

pub struct PairedChecker {
    pub read_length: usize,
    pub min_quality: u8,
    pub min_dist: i32,
    pub max_dist: i32,
    pub force_paired: bool,
}

impl RecordCheck for PairedChecker {
    fn valid(&self, record: &Record) -> bool {
        // check single read conditions
        if record.is_unmapped() || record.seq().len() == self.read_length || record.mapq() >= self.min_quality {
            return false;
        }
        // mandatory paired condition
        if ( !record.is_paired() || record.is_mate_unmapped() || record.tid() != record.mtid() ) && self.force_paired {
            return false;
        }
        // check pair distance 
        if record.is_paired() {
            let dist = (record.pos() - record.mpos()).abs() + self.read_length as i32;
            if dist < self.min_dist || dist > self.max_dist {
                return false;
            }
        }
        true
    }
}
