use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::info;
use paraseq::{
    fastq,
    fastx::Record,
    parallel::{ParallelProcessor, ParallelReader, ProcessError},
};
use parking_lot::Mutex;
// use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

type ArcHash = Arc<Mutex<HashMap<usize, usize>>>;

#[derive(Clone)]
struct Length {
    len: HashMap<usize, usize>,
    total: ArcHash,
}

impl Length {
    pub fn new() -> Self {
        Self {
            len: HashMap::new(),
            total: Arc::new(Mutex::new(HashMap::new())),
        }
    }
}

impl ParallelProcessor for Length {
    fn process_record<Rf: Record>(&mut self, record: Rf) -> Result<(), ProcessError> {
        let rlen = record.seq().len();
        *self.len.entry(rlen).or_insert(0usize) += 1;

        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<(), ProcessError> {
        let mut total = self.total.lock();
        for (k, v) in self.len.iter() {
            *total.entry(*k).or_insert(0usize) += *v;
        }

        // reset the len HashMap for the next batch
        self.len.clear();
        Ok(())
    }
}

pub fn fq_length(
    file: Option<&String>,
    rev: bool,
    ncpu: usize,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    let length = Length::new();
    fq_reader.process_parallel(length.clone(), ncpu)?;

    let mut fo = file_writer(out, compression_level, stdout_type)?;

    let reads_len = length.total.lock();
    let mut sort_len: Vec<(&usize, &usize)> = reads_len.iter().collect();

    if rev {
        // sort_len.par_sort_by(|x, y| y.0.cmp(x.0));
        sort_len.sort_by(|x, y| y.0.cmp(x.0));
    } else {
        // sort_len.par_sort_by_key(|x| x.0);
        sort_len.sort_by_key(|x| x.0);
    }

    fo.write_all("lenth\tcount\n".as_bytes())?;
    for (k, v) in sort_len.iter() {
        fo.write_all(format!("{}\t{}\n", k, v).as_bytes())?;
    }
    fo.flush()?;
    info!(
        "total scan reads number: {}",
        reads_len.iter().map(|(_, v)| *v).sum::<usize>()
    );

    Ok(())
}
