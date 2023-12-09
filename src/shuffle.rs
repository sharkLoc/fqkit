use crate::utils::*;
use bio::io::fastq;
use std::io::Result;
use log::*;
use std::time::Instant;
use rand::prelude::*;
use rand_pcg::Pcg64;


pub fn shuffle_fastq(
    file: &Option<&str>,
    seed: u64,
    out: &Option<&str>,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        if let Some(file) = file {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
        info!("rand seed: {}",seed);
    }
    let start = Instant::now();
    let mut rng = Pcg64::seed_from_u64(seed);
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    
    let mut vec_reads = vec![];
    for rec in fq_reader.records().flatten() {
        vec_reads.push(rec);
    }
    
    if !quiet { info!("all records has been readed into memory, start shuffle ..."); }
    vec_reads.shuffle(&mut rng);
    if !quiet { info!("shuffle done, start write to output ..."); }

    let mut fq_writer = file_writer(out).map(fastq::Writer::new)?;
    for rec in vec_reads {
        fq_writer.write_record(&rec)?;
    }
    fq_writer.flush()?;

    if !quiet {
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}