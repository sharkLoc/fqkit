use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;

pub fn interleaved(
    file1: &String,
    file2: &String,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {

    let mut num = 0usize;
    let fq1_reader = fastq::Reader::new(file_reader(Some(file1))?);
    let fq2_reader = fastq::Reader::new(file_reader(Some(file2))?);
    info!("reading from file: {}", file1);
    info!("reading from file: {}", file2);

    let mut fq_writer = fastq::Writer::new(file_writer(out, compression_level, stdout_type)?);
    for (rec1, rec2) in fq1_reader
        .records()
        .map_while(Result::ok)
        .zip(fq2_reader.records().map_while(Result::ok))
    {
        num += 1;
        fq_writer.write(rec1.id(), rec1.desc(), rec1.seq(), rec1.qual())?;
        fq_writer.write(rec2.id(), rec2.desc(), rec2.seq(), rec2.qual())?;
    }
    fq_writer.flush()?;

    info!("total PE reads number: {}", num);
    Ok(())
}
