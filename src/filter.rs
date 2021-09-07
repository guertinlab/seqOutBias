use htslib::bam::record::Record;

pub trait RecordCheck {
    fn valid(&self, rec: &Record) -> bool;
    fn vir_pos(&self, rec: &Record) -> i32;
}

pub struct SingleChecker {
    pub tail_edge: bool,
    pub exact_length: bool,
    pub read_length: usize,
    pub min_quality: u8,
}

fn vir_pos_common(read_length: usize, tail_edge: bool, exact_length: bool, rec: &Record) -> i32 {
    if tail_edge {
        // instead of right edge of read (5') report left edge (3')
        if rec.is_reverse() {
            rec.pos() - (rec.seq().len() as i32) + 1
        } else {
            rec.pos() + (rec.seq().len() as i32) - 1
        }
    } else {
        if !exact_length && rec.is_reverse() {
            rec.pos() + (rec.seq().len() as i32) - (read_length as i32)
        } else {
            rec.pos()
        }
    }
}

impl RecordCheck for SingleChecker {
    fn valid(&self, record: &Record) -> bool {
        !record.is_unmapped() && (!self.exact_length || record.seq().len() == self.read_length) && record.mapq() >= self.min_quality
    }

    fn vir_pos(&self, rec: &Record) -> i32 {
        vir_pos_common(self.read_length, self.tail_edge, self.exact_length, rec)
    }
}

#[derive(Copy, Clone)]
pub enum PairPosition {
    First,
    Last
}

pub struct PairedChecker {
    pub tail_edge: bool,
    pub exact_length: bool,
    pub read_length: usize,
    pub min_quality: u8,
    pub min_dist: i32,
    pub max_dist: i32,
    pub force_paired: bool,
    pub max_distance: bool,
    pub select_pair: Option<PairPosition>
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
        // filter for specific pair in paired reads
        match self.select_pair {
            Some( side ) => {
                match side {
                    PairPosition::First => {
                        if !record.is_first_in_template() {
                            return false;
                        }
                    }
                    PairPosition::Last => {
                        if !record.is_last_in_template() {
                            return false;
                        }
                    }
                }
            }
            _ => {}
        }
        true
    }

    fn vir_pos(&self, rec: &Record) -> i32 {
        vir_pos_common(self.read_length, self.tail_edge, self.exact_length, rec)
    }
}
