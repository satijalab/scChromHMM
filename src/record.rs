use std::collections::HashMap;
use std::ops::Range;

use crate::config::{ProbT, RangeT};

#[derive(Debug, Clone)]
pub struct Record<R> {
    range: Range<RangeT>,
    id: R,
}

impl<R> Record<R>
where
    R: Copy,
{
    pub fn id(&self) -> R {
        self.id
    }

    pub fn range(&self) -> &Range<RangeT> {
        &self.range
    }

    pub fn new_with_id(range: &Range<RangeT>, id: R) -> Record<R> {
        Record {
            range: Range {
                start: range.start,
                end: range.end,
            },
            id,
        }
    }
}

impl Record<u64> {
    pub fn from_string(
        data: String,
        assay_cells: &HashMap<u64, HashMap<u32, ProbT>>,
    ) -> Option<Record<u64>> {
        let (mut start, mut end) = (0, 0);
        let mut cb: &str = "chr0";
        let mut id: u8 = 0;
        for (index, text) in data.split_whitespace().enumerate() {
            match index {
                1 => start = text.parse::<u32>().unwrap(),
                2 => end = text.parse::<u32>().unwrap(),
                3 => {
                    let toks: Vec<&str> = text.split('-').collect();
                    cb = toks[0];
                    id = toks
                        .get(1)
                        .unwrap_or_else(|| panic!("{:?}", toks.get(0)))
                        .parse::<u8>()
                        .unwrap();
                }
                4 => break,
                _ => (),
            }
        }

        let id = match id {
            1 => carina::barcode::cb_string_to_u64_with_id(cb.as_bytes(), 16, 1).unwrap(),
            2 => carina::barcode::cb_string_to_u64_with_id(cb.as_bytes(), 16, 2).unwrap(),
            _ => unreachable!(),
        };

        match assay_cells.contains_key(&id) {
            true => {
                // 1-offset
                let range = Range {
                    start: start - 1,
                    end: end,
                };
                Some(Record { range, id })
            }
            false => None,
        }
    }
}

////////////////////////////////////////////
/// Cell Records
////////////////////////////////////////////
#[derive(Debug)]
pub struct CellRecords<R> {
    records: Vec<Record<R>>,
}

impl CellRecords<ProbT> {
    pub fn new(records: Vec<Record<ProbT>>) -> CellRecords<ProbT> {
        CellRecords { records }
    }

    pub fn records(&self) -> &Vec<Record<ProbT>> {
        &self.records
    }
}

////////////////////////////////////////////
/// Assay Records
////////////////////////////////////////////
#[derive(Debug)]
pub struct AssayRecords<R> {
    records: Vec<CellRecords<R>>,
}

impl AssayRecords<ProbT> {
    pub fn new(records: Vec<CellRecords<ProbT>>) -> AssayRecords<ProbT> {
        AssayRecords { records }
    }

    pub fn get_cell_records(&self, cell_id: usize) -> Option<&CellRecords<ProbT>> {
        self.records.get(cell_id)
    }
}

////////////////////////////////////////////
/// Experiment
////////////////////////////////////////////
#[derive(Debug)]
pub struct Experiment<ProbT> {
    records: Vec<AssayRecords<ProbT>>,
}

impl Experiment<ProbT> {
    pub fn new(records: Vec<AssayRecords<ProbT>>) -> Experiment<ProbT> {
        Experiment { records }
    }

    pub fn get_cell_data(&self, cell_id: usize) -> Vec<&CellRecords<ProbT>> {
        self.records
            .iter()
            .map(|x| x.get_cell_records(cell_id).unwrap())
            .collect()
    }
}
