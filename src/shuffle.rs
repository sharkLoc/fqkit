use crate::utils::*;
use bio::io::fastq;
use anyhow::Result;
use log::*;
use std::time::Instant;
use rand::prelude::*;
use rand_pcg::Pcg64;


pub fn shuffle_fastq(
    file: &Option<&str>,
    seed: u64,
    out: &Option<&str>,
    compression_level: u32,
) -> Result<()> {
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("rand seed: {}",seed);
    let start = Instant::now();

    let mut rng = Pcg64::seed_from_u64(seed);
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    
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

    info!("time elapsed is: {:?}",start.elapsed());
    Ok(())
}
