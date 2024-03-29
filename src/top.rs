use anyhow::Result;
use crate::utils::*;
use bio::io::fastq;
use std::time::Instant;
use log::*;

pub fn top_n_records(
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
    info!("get top {} records",number);
    let start = Instant::now();

    let fp = fastq::Reader::new(file_reader(input)?);
    let mut fo = fastq::Writer::new(file_writer(output, compression_level)?);
    for rec in fp.records().take(number).flatten() {
        fo.write_record(&rec)?;
    }
    fo.flush()?;

    info!("time elapsed is: {:?}",start.elapsed());
    Ok(())
}
