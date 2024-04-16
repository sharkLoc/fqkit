use crate::utils::*;
use anyhow::Result;
use bio::io::fastq;
use log::*;
use rand::prelude::*;
use rand_pcg::Pcg64;
use std::time::Instant;

pub fn shuffle_fastq(
    file: Option<&String>,
    seed: u64,
    out: Option<&String>,
    compression_level: u32,
) -> Result<()> {
    let start = Instant::now();

    let mut rng = Pcg64::seed_from_u64(seed);
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("rand seed: {}", seed);
    
    let mut vec_reads = vec![];
    for rec in fq_reader.records().flatten() {
        vec_reads.push(rec);
    }

    info!("all records has been readed into memory, start shuffle ...");
    vec_reads.shuffle(&mut rng);
    info!("shuffle done, start write to output ...");

    let mut fq_writer = file_writer(out, compression_level).map(fastq::Writer::new)?;
    for rec in vec_reads {
        fq_writer.write_record(&rec)?;
    }
    fq_writer.flush()?;

    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
