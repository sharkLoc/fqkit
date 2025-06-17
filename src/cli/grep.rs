use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{error, info};
use paraseq::{
    fastq,
    fastx::Record,
    parallel::{ParallelProcessor, ParallelReader, ProcessError},
};
use parking_lot::Mutex;
use std::io::{BufRead, Write};
use std::sync::Arc;

#[derive(Clone)]
struct Grepid {
    writer: Arc<Mutex<Box<dyn Write + Send>>>,
    count: usize,
    clount_all: Arc<Mutex<usize>>,
    full_name: bool,
    list: Vec<Vec<u8>>,
    buffer: Vec<u8>,
}

impl Grepid {
    pub fn new(
        writer: Box<dyn Write + Send>,
        count: usize,
        clount_all: usize,
        full_name: bool,
        list: Vec<Vec<u8>>,
        buffer: Vec<u8>,
    ) -> Self {
        Self {
            writer: Arc::new(Mutex::new(writer)),
            count,
            clount_all: Arc::new(Mutex::new(clount_all)),
            full_name,
            list,
            buffer,
        }
    }

    pub fn write_record<Rf: Record>(&mut self, record: Rf) -> std::io::Result<()> {
        self.count += 1;
        write_record(
            &mut self.buffer,
            record.id(),
            record.seq(),
            record.qual().unwrap(),
        )?;
        Ok(())
    }
}

impl ParallelProcessor for Grepid {
    fn process_record<Rf: Record>(&mut self, record: Rf) -> Result<(), ProcessError> {
        let name = if self.full_name {
            record.id().to_vec()
        } else {
            record
                .id_str()
                .split_whitespace()
                .next()
                .unwrap()
                .as_bytes()
                .to_vec()
        };

        if self.list.contains(&name) {
            self.write_record(record)?;
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        let mut writer = self.writer.lock();

        writer.write_all(&self.buffer)?;
        writer.flush()?;

        let mut clount_all = self.clount_all.lock();
        *clount_all += self.count;

        // reset for next batch
        self.count = 0;
        self.buffer.clear();
        Ok(())
    }
}

pub fn grep_fastq(
    fq: Option<&String>,
    list: &String,
    full_name: bool,
    ncpu: usize,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut ids = vec![];
    let fp_id = file_reader(Some(list))?;
    info!("reading reads id from file: {}", list);
    for id in fp_id.lines().map_while(Result::ok) {
        ids.push(id.as_bytes().to_vec());
    }
    if ids.is_empty() {
        error!("{}", FqkitError::EmptyFile(list.to_string()));
        std::process::exit(1);
    }

    let fq_writer = file_writer(out, compression_level, stdout_type)?;

    let grepid = Grepid::new(fq_writer, 0, 0, full_name, ids, vec![]);

    let fq_reader = file_reader(fq).map(fastq::Reader::new)?;
    fq_reader.process_parallel(grepid.clone(), ncpu)?;

    if let Some(out) = out {
        info!("reads write to file: {}", out);
    } else {
        info!("reads write to stdout");
    }

    info!("total reads matched number: {}", *grepid.clount_all.lock());
    Ok(())
}
