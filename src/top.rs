use std::io::Result;
use crate::utils::*;
use bio::io::fastq;
use std::time::Instant;
use log::*;

pub fn top_n_records(
    input: &Option<&str>,
    number: usize,
    output: &Option<&str>,
) -> Result<()> {
    info!("reading from file: {}",input.unwrap());
    info!("get top {} records",number);
    let start = Instant::now();

    let fp = fastq::Reader::new(file_reader(input)?);
    let mut fo = fastq::Writer::new(file_writer(output)?);
    for rec in fp.records().take(number).flatten() {
        fo.write_record(&rec)?;
    }
  
    info!("time elapsed is: {:?}",start.elapsed());
    Ok(())
}