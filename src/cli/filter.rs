use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{error, info};
use paraseq::{
    fastq,
    fastx::Record,
    parallel::{PairedParallelProcessor, PairedParallelReader, ProcessError},
};
use parking_lot::Mutex;
use std::{io::Write, sync::Arc};

#[derive(Clone)]
struct FilterSeq {
    count_n: usize,
    length: usize,
    complexity: usize,
    average_qual: u8,
    phred: u8,
    buffer1: Vec<u8>,
    buffer2: Vec<u8>,
    failed_buffer: Vec<u8>,
    pe_ok: usize,
    pe_fail: usize,
    total_pe_ok: Arc<Mutex<usize>>,
    total_pe_fail: Arc<Mutex<usize>>,
    writer1: Arc<Mutex<Box<dyn Write + Send>>>,
    writer2: Arc<Mutex<Box<dyn Write + Send>>>,
    failed_writer: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
}

impl FilterSeq {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        count_n: usize,
        length: usize,
        complexity: usize,
        average_qual: u8,
        phred: u8,
        buffer1: Vec<u8>,
        buffer2: Vec<u8>,
        failed_buffer: Vec<u8>,
        pe_ok: usize,
        pe_fail: usize,
        output1: Box<dyn Write + Send>,
        output2: Box<dyn Write + Send>,
        failed_out: Option<Box<dyn Write + Send>>,
    ) -> Self {
        Self {
            count_n,
            length,
            complexity,
            average_qual,
            phred,
            buffer1,
            buffer2,
            failed_buffer,
            pe_ok,
            pe_fail,
            total_pe_ok: Arc::new(Mutex::new(0)),
            total_pe_fail: Arc::new(Mutex::new(0)),
            writer1: Arc::new(Mutex::new(output1)),
            writer2: Arc::new(Mutex::new(output2)),
            failed_writer: failed_out.map(|failed_out| Arc::new(Mutex::new(failed_out))),
        }
    }

    pub fn qc_base_n(&self, seq: &[u8]) -> bool {
        seq.iter().filter(|v| *v == &b'N').count() <= self.count_n
    }

    pub fn qc_phread_mean(&self, qual: &[u8]) -> bool {
        phred_mean(qual, self.phred) >= self.average_qual
    }

    pub fn qc_length(&self, seq: &[u8]) -> bool {
        seq.len() >= self.length
    }

    pub fn qc_complx(&self, seq: &[u8]) -> bool {
        let complex = seq
            .iter()
            .skip(1)
            .zip(seq.iter())
            .filter(|(q1, q2)| q1 != q2)
            .count() as f64
            / (seq.len() - 1) as f64;
        complex as usize >= self.complexity
    }

    pub fn qc_all(&self, seq: &[u8], qual: &[u8]) -> bool {
        self.qc_base_n(seq)
            && self.qc_phread_mean(qual)
            && self.qc_complx(seq)
            && self.qc_length(seq)
    }

    pub fn write_record1<Rf: Record>(&mut self, record: Rf) -> std::io::Result<()> {
        write_record(
            &mut self.buffer1,
            record.id(),
            record.seq(),
            record.qual().unwrap(),
        )?;
        // self.buffer1.write_all(b"@")?;
        // self.buffer1.extend_from_slice(record.id());
        // self.buffer1.write_all(b"\n")?;
        // self.buffer1.extend_from_slice(record.seq());
        // self.buffer1.write_all(b"\n+\n")?;
        // self.buffer1.extend_from_slice(record.qual().unwrap());
        // self.buffer1.write_all(b"\n")?;
        Ok(())
    }

    pub fn write_record2<Rf: Record>(&mut self, record: Rf) -> std::io::Result<()> {
        write_record(
            &mut self.buffer2,
            record.id(),
            record.seq(),
            record.qual().unwrap(),
        )?;
        // self.buffer2.write_all(b"@")?;
        // self.buffer2.extend_from_slice(record.id());
        // self.buffer2.write_all(b"\n")?;
        // self.buffer2.extend_from_slice(record.seq());
        // self.buffer2.write_all(b"\n+\n")?;
        // self.buffer2.extend_from_slice(record.qual().unwrap());
        // self.buffer2.write_all(b"\n")?;
        Ok(())
    }

    pub fn write_record_fail<Rf: Record>(&mut self, record: Rf) -> std::io::Result<()> {
        write_record(
            &mut self.failed_buffer,
            record.id(),
            record.seq(),
            record.qual().unwrap(),
        )?;
        // self.failed_buffer.write_all(b"@")?;
        // self.failed_buffer.extend_from_slice(record.id());
        // self.failed_buffer.write_all(b"\n")?;
        // self.failed_buffer.extend_from_slice(record.seq());
        // self.failed_buffer.write_all(b"\n+\n")?;
        // self.failed_buffer.extend_from_slice(record.qual().unwrap());
        // self.failed_buffer.write_all(b"\n")?;
        Ok(())
    }
}

impl PairedParallelProcessor for FilterSeq {
    fn process_record_pair<Rf: Record>(&mut self, rec1: Rf, rec2: Rf) -> Result<(), ProcessError> {
        if self.qc_all(rec1.seq(), rec1.qual().unwrap())
            && self.qc_all(rec2.seq(), rec2.qual().unwrap())
        {
            self.write_record1(rec1)?;
            self.write_record2(rec2)?;
            self.pe_ok += 1;
        } else {
            if self.failed_writer.is_some() {
                self.write_record_fail(rec1)?;
                self.write_record_fail(rec2)?;
            }
            self.pe_fail += 1;
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<(), ProcessError> {
        let mut writer1 = self.writer1.lock();
        let mut writer2 = self.writer2.lock();
        let mut pe_ok = self.total_pe_ok.lock();
        let mut pe_fail = self.total_pe_fail.lock();

        *pe_ok += self.pe_ok;
        *pe_fail += self.pe_fail;

        writer1.write_all(&self.buffer1)?;
        writer2.write_all(&self.buffer2)?;
        writer1.flush()?;
        writer2.flush()?;

        if let Some(failed_writer) = &self.failed_writer {
            let mut writer = failed_writer.lock();
            writer.write_all(&self.failed_buffer)?;
            writer.flush()?;
            // reset for next batch
            self.failed_buffer.clear();
        }

        // reset for next batch
        self.buffer1.clear();
        self.buffer2.clear();
        self.pe_ok = 0;
        self.pe_fail = 0;
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub fn filter_fastq(
    read1: &String,
    read2: &String,
    nbase: usize,
    length: usize,
    complexity: u32,
    average_qual: u8,
    phred: u8,
    ncpu: usize,
    failed: Option<&String>,
    out1: &String,
    out2: &String,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    if ![33u8, 64u8].contains(&phred) {
        error!("{}", FqkitError::InvalidPhredValue);
        std::process::exit(1);
    }

    let fq_reader1 = file_reader(Some(read1)).map(fastq::Reader::new)?;
    let fq_reader2 = file_reader(Some(read2)).map(fastq::Reader::new)?;
    let out_writer1 = file_writer(Some(out1), compression_level, stdout_type)?;
    let out_writer2 = file_writer(Some(out2), compression_level, stdout_type)?;

    let failed_writer = if let Some(failed) = failed {
        Some(file_writer(Some(failed), compression_level, stdout_type)?)
    } else {
        None
    };

    let complex = complexity as usize;
    let (pe_ok, pe_fail) = (0usize, 0usize);
    let buffer1 = vec![];
    let buffer2 = vec![];
    let failed_buffer = vec![];

    let filters = FilterSeq::new(
        nbase,
        length,
        complex,
        average_qual,
        phred,
        buffer1,
        buffer2,
        failed_buffer,
        pe_ok,
        pe_fail,
        out_writer1,
        out_writer2,
        failed_writer,
    );
    // run the filter
    fq_reader1.process_parallel_paired(fq_reader2, filters.clone(), ncpu)?;

    let pe_ok = filters.total_pe_ok.lock();
    let pe_fail = filters.total_pe_fail.lock();
    info!("total clean pe reads number (r1+r2): {}", *pe_ok * 2);
    info!("total failed pe reads number (r1+r2): {}", *pe_fail * 2);
    Ok(())
}

fn phred_mean(qual: &[u8], phred: u8) -> u8 {
    let ave_error = qual
        .iter()
        .map(|x| 10.0f64.powf((x - phred) as f64 / -10.0))
        .sum::<f64>()
        / qual.len() as f64;

    (-10.0f64 * ave_error.log10()).round() as u8
}
