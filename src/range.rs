use std::io::Result;
use crate::utils::*;
use bio::io::fastq;
use std::time::Instant;
use log::*;

pub fn range_fastq(
    input: &Option<&str>,
    skip: usize,
    take: usize,
    output: &Option<&str>,
    compression_level: u32,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        if let Some(file) = input {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
        info!("skip first {} records", skip);
        info!("get {} records", take);
    }
    let start = Instant::now();

    let fp_reader = file_reader(input).map(fastq::Reader::new)?;
    let mut fp_writer = file_writer(output, compression_level).map(fastq::Writer::new)?;
    
    for rec in fp_reader.records().skip(skip).take(take).flatten() {
        fp_writer.write_record(&rec)?;
    }
    fp_writer.flush()?;

    if !quiet {
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}