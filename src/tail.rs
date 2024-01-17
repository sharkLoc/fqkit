use anyhow::Result;
use crate::utils::*;
use bio::io::fastq;
use std::time::Instant;
use log::*;

pub fn tail_n_records(
    input: &Option<&str>,
    number: usize,
    output: &Option<&str>,
    compression_level: u32,
) -> Result<()> {
    if let Some(file) = input {
        info!("reading from file: {}",file);
    } else {
        info!("reading from stdin");
    }
    info!("get tail {} records",number);
    let start = Instant::now();

    let fp = fastq::Reader::new(file_reader(input)?);
    let mut fo = fastq::Writer::new(file_writer(output, compression_level)?);
    let mut total = 0usize;

    for _ in fp.records() { total += 1; }
    info!("fastq file total reads number: {}",total);
    let skip_n = total - number;
    
    let fp2 = fastq::Reader::new(file_reader(input)?);
    for rec in fp2.records().skip(skip_n).flatten() {
        fo.write_record(&rec)?;
    }
    fo.flush()?;

    info!("time elapsed is: {:?}",start.elapsed());
    Ok(())
}