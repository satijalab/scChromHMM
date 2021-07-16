use crate::config::ProbT;
use crate::record::{CellRecords, Record};
use std::ops::Range;

use rust_htslib::tbx::{self, Read};
use std::collections::HashMap;
use std::path::PathBuf;

pub struct Fragment {
    _filepath: PathBuf,
    reader: tbx::Reader,
}

impl Fragment {
    pub fn from_pathbuf(filepath: PathBuf) -> Fragment {
        let tbx_reader = tbx::Reader::from_path(&filepath)
            .unwrap_or_else(|_| panic!("Could not open {:?}", filepath));

        Fragment {
            _filepath: filepath,
            reader: tbx_reader,
        }
    }

    pub fn _filepath(&self) -> &str {
        self._filepath.to_str().unwrap()
    }

    pub fn tid(&self, seqname: &str) -> u64 {
        match self.reader.tid(seqname) {
            Ok(tid) => tid,
            Err(_) => panic!("Could not resolve to contig ID"),
        }
    }

    pub fn fetch(
        &mut self,
        tid: u64,
        region: &Range<u32>,
        assay_cells: &HashMap<u64, HashMap<u32, ProbT>>,
        num_common_cells: usize,
    ) -> Vec<CellRecords<ProbT>> {
        // Set region to fetch.
        self.reader
            .fetch(tid, region.start as u64, region.end as u64)
            .expect("Could not seek to fetch region");

        let all_records: Vec<Record<u64>> = self
            .reader
            .records()
            .map(|x| String::from_utf8(x.unwrap()).expect("UTF8 conversion error"))
            .map(|x| Record::from_string(x, assay_cells))
            .flatten()
            .collect();

        let mut norm: HashMap<u32, u32> = HashMap::new();
        for (_, val) in assay_cells {
            for (cb, _) in val {
                let ct = norm.entry(*cb).or_insert(0);
                *ct += 1;
            }
        }

        let mut cell_records: Vec<Vec<Record<ProbT>>> = vec![Vec::new(); num_common_cells];
        for record in all_records {
            let cb = record.id();
            if let Some(dict) = assay_cells.get(&cb) {
                for (&cell_id, &prob) in dict {
                    let nprob = prob / *norm.get(&cell_id).unwrap() as ProbT;
                    let new_record = Record::new_with_id(&record.range(), nprob);
                    cell_records[cell_id as usize].push(new_record);
                }
            };
        }

        let cell_records: Vec<CellRecords<ProbT>> =
            cell_records.into_iter().map(CellRecords::new).collect();

        cell_records
    }
}
