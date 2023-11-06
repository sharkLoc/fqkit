use crate::utils::*;
use bio::io::fastq;
use anyhow::{Error, Ok};
use log::*;
use std::time::Instant;


pub fn interleaved(
    file1: &Option<&str>,
    file2: &Option<&str>,
    out: &Option<&str>,
    quiet: bool,
) -> Result<(), Error> {
    
    if !quiet {
        info!("reading from file: {}", file1.unwrap());
        info!("reading from file: {}", file2.unwrap());
    }
    let start = Instant::now();

    let mut num = 0usize;
    let fq1_reader = fastq::Reader::new(file_reader(file1)?);
    let fq2_reader = fastq::Reader::new(file_reader(file2)?);
    let mut fq_writer = fastq::Writer::new(file_writer(out)?);
    
    for (rec1, rec2) in fq1_reader.records().flatten().zip(fq2_reader.records().flatten()) {
        num += 1;
        fq_writer.write(rec1.id(), rec1.desc(), rec1.seq(), rec1.qual())?;
        fq_writer.write(rec2.id(), rec2.desc(), rec2.seq(), rec2.qual())?;
    }

    if !quiet {
        info!("total PE reads number: {}",num);
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}
