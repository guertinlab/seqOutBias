use htslib::bam::record::Record;

pub trait RecordCheck {
    fn valid(&self, rec: &Record) -> bool;
    fn vir_pos(&self, rec: &Record) -> u32;
}

pub struct SingleChecker {
    pub exact_length: bool,
    pub read_length: usize,
    pub min_quality: u8,
}

impl RecordCheck for SingleChecker {
    fn valid(&self, record: &Record) -> bool {
        !record.is_unmapped() && (!self.exact_length || record.seq().len() == self.read_length) && record.mapq() >= self.min_quality
    }

    fn vir_pos(&self, rec: &Record) -> u32 {
        if self.exact_length && rec.is_reverse() {
            let aux = rec.pos() + (rec.seq().len() - self.read_length) as i32;
            if aux >= 0 {
                aux as u32
            }
            else {
                0
            }
        }
        else {
            rec.pos() as u32
        }
    }
}

pub struct PairedChecker {
    pub exact_length: bool,
    pub read_length: usize,
    pub min_quality: u8,
    pub min_dist: i32,
    pub max_dist: i32,
    pub force_paired: bool,
    pub max_distance: bool,
}

impl RecordCheck for PairedChecker {
    fn valid(&self, record: &Record) -> bool {
        // check single read conditions
        if record.is_unmapped() || (self.exact_length && record.seq().len() != self.read_length) || record.mapq() < self.min_quality {
            return false;
        }
        // mandatory paired condition
        if ( !record.is_paired() || record.is_mate_unmapped() || record.tid() != record.mtid() ) && self.force_paired {
            return false;
        }
        // check pair distance 
        if record.is_paired() && self.max_distance {
            let dist = (record.pos() - record.mpos()).abs() + self.read_length as i32; // TODO FIX!
            if dist < self.min_dist || dist > self.max_dist {
                return false;
            }
        }
        true
    }

    fn vir_pos(&self, rec: &Record) -> u32 {
        if self.exact_length && rec.is_reverse() {
            let aux = rec.pos() + (rec.seq().len() - self.read_length) as i32;
            if aux >= 0 {
                aux as u32
            }
            else {
                0
            }
        }
        else {
            rec.pos() as u32
        }
    }
}
