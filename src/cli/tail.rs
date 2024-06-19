use crate::utils::*;
use anyhow::Result;
use bio::io::fastq;
use log::*;
use std::time::Instant;

pub fn tail_n_records(
    input: Option<&String>,
    number: usize,
    rdc: bool,
    output: Option<&String>,
    compression_level: u32,
) -> Result<()> {
    let start = Instant::now();

    let fp = fastq::Reader::new(file_reader(input)?);
    if let Some(file) = input {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("get tail {} records", number);

    let mut fo = fastq::Writer::new(file_writer(output, compression_level)?);
    if rdc {
        let mut total = 0usize;
        for _ in fp.records() {
            total += 1;
        }
        info!("fastq file total reads number: {}", total);

        let skip_n = total - number;
        let fp2 = fastq::Reader::new(file_reader(input)?);
        for rec in fp2.records().skip(skip_n).map_while(Result::ok) {
            fo.write_record(&rec)?;
        }
    } else {
        let mut tail = vec![];
        for rec in fp.records().map_while(Result::ok) {
            tail.push(rec);
        }
        for rec in tail.iter().rev().take(number).rev() {
            fo.write_record(rec)?;
        }
    }
    fo.flush()?;

    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
