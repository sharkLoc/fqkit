use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::*;
use paraseq::{
    fastq,
    fastx::Record,
    parallel::{ParallelProcessor, ParallelReader, ProcessError},
};
use parking_lot::Mutex;
use regex::RegexBuilder;
use std::io::Write;
use std::sync::Arc;

#[derive(Clone)]
struct Search {
    pat: regex::Regex,
    invert_match: bool,
    count: usize,
    buffer: Vec<u8>,
    tatal_count: Arc<Mutex<usize>>,
    writer: Arc<Mutex<Box<dyn Write + Send>>>,
}

impl Search {
    pub fn new(pat: &str, case: bool, invert_match: bool, writer: Box<dyn Write + Send>) -> Self {
        let re = RegexBuilder::new(pat)
            .case_insensitive(case)
            .build()
            .unwrap();
        Search {
            pat: re,
            invert_match,
            count: 0,
            buffer: Vec::new(),
            tatal_count: Arc::new(Mutex::new(0)),
            writer: Arc::new(Mutex::new(writer)),
        }
    }
}

impl ParallelProcessor for Search {
    fn process_record<Rf: Record>(&mut self, record: Rf) -> Result<(), ProcessError> {
        let fq_str = record.seq_str();
        if (self.invert_match && !self.pat.is_match(fq_str))
            || (!self.invert_match && self.pat.is_match(fq_str))
        {
            self.count += 1;
            write_record(
                &mut self.buffer,
                record.id(),
                record.seq(),
                record.qual().unwrap(),
            )?;
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<(), ProcessError> {
        let mut writer = self.writer.lock();
        writer.write_all(&self.buffer)?;

        let mut total_count = self.tatal_count.lock();
        *total_count += self.count;

        self.count = 0;
        self.buffer.clear();
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub fn search_fq(
    fq: Option<&String>,
    pat: &str,
    case: bool,
    invert_match: bool,
    out: Option<&String>,
    ncpu: usize,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let fq_reader = file_reader(fq).map(fastq::Reader::new).unwrap();
    info!("regex pattern is: {}", pat);

    let fo = file_writer(out, compression_level, stdout_type)?;
    if let Some(out) = out {
        info!("reads write to file: {}", out);
    } else {
        info!("reads write to stdout");
    }

    let processor = Search::new(pat, case, invert_match, fo);
    fq_reader.process_parallel(processor.clone(), ncpu)?;

    info!("total reads number: {}", processor.tatal_count.lock());
    Ok(())
}
