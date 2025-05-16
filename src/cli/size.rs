use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use paraseq::{
    fastq,
    fastx::Record,
    parallel::{ParallelProcessor, ParallelReader, ProcessError},
};
use parking_lot::Mutex;
use std::sync::Arc;

#[derive(Clone, Default)]
struct Base {
    a: usize,
    t: usize,
    g: usize,
    c: usize,
    n: usize,
    read: usize,
    total_a: Arc<Mutex<usize>>,
    total_t: Arc<Mutex<usize>>,
    total_g: Arc<Mutex<usize>>,
    total_c: Arc<Mutex<usize>>,
    total_n: Arc<Mutex<usize>>,
    total_bases: Arc<Mutex<usize>>,
    total_reads: Arc<Mutex<usize>>,
}

impl Base {
    pub fn get_total_a(&self) -> usize {
        *self.total_a.lock()
    }

    pub fn get_total_t(&self) -> usize {
        *self.total_t.lock()
    }

    pub fn get_total_g(&self) -> usize {
        *self.total_g.lock()
    }

    pub fn get_total_c(&self) -> usize {
        *self.total_c.lock()
    }

    pub fn get_total_n(&self) -> usize {
        *self.total_n.lock()
    }

    pub fn get_total_bases(&self) -> usize {
        *self.total_bases.lock()
    }

    pub fn get_total_reads(&self) -> usize {
        *self.total_reads.lock()
    }
}

impl ParallelProcessor for Base {
    fn process_record<Rf: Record>(&mut self, record: Rf) -> Result<(), ProcessError> {
        for nt in record.seq().iter() {
            match *nt {
                b'A' => self.a += 1,
                b'T' => self.t += 1,
                b'G' => self.g += 1,
                b'C' => self.c += 1,
                b'N' => self.n += 1,
                _ => unreachable!(),
            }
        }
        self.read += 1;

        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<(), ProcessError> {
        let mut total_a = self.total_a.lock();
        let mut total_t = self.total_t.lock();
        let mut total_g = self.total_g.lock();
        let mut total_c = self.total_c.lock();
        let mut total_n = self.total_n.lock();
        let mut total_bases = self.total_bases.lock();
        let mut total_reads = self.total_reads.lock();

        *total_a += self.a;
        *total_t += self.t;
        *total_g += self.g;
        *total_c += self.c;
        *total_n += self.n;
        *total_bases += self.a + self.t + self.g + self.c + self.n;
        *total_reads += self.read;

        // Reset the counts for the next batch
        self.a = 0;
        self.t = 0;
        self.g = 0;
        self.c = 0;
        self.n = 0;
        self.read = 0;

        Ok(())
    }
}

pub fn size_fastq(
    fq: Option<&String>,
    ncpu: usize,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let processor = Base::default();

    let fq_reader = fastq::Reader::new(file_reader(fq)?);
    fq_reader.process_parallel(processor.clone(), ncpu)?;

    let mut fo = file_writer(out, compression_level, stdout_type)?;
    fo.write_all(
        format!(
            "reads:{}\tbases:{}\tA:{}\tT:{}\tG:{}\tC:{}\tN:{}\n",
            processor.get_total_reads(),
            processor.get_total_bases(),
            processor.get_total_a(),
            processor.get_total_t(),
            processor.get_total_g(),
            processor.get_total_c(),
            processor.get_total_n()
        )
        .as_bytes(),
    )?;

    Ok(())
}
