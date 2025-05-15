use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::info;
use rand::prelude::*;
use rand_pcg::Pcg64;

pub fn shuffle_fastq(
    file: Option<&String>,
    seed: u64,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut rng = Pcg64::seed_from_u64(seed);
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    info!("rand seed: {}", seed);

    let mut vec_reads = vec![];
    for rec in fq_reader.records().map_while(Result::ok) {
        vec_reads.push(rec);
    }

    info!("all records has been readed into memory, start shuffle ...");
    vec_reads.shuffle(&mut rng);
    info!("shuffle done, start write to output ...");

    let mut fq_writer = file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;
    for rec in vec_reads {
        fq_writer.write_record(&rec)?;
    }
    fq_writer.flush()?;

    Ok(())
}
