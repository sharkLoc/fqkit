use crate::utils::*;
use anyhow::Result;
use bio::io::fastq;
use log::*;

pub fn range_fastq(
    input: Option<&String>,
    skip: usize,
    take: usize,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<()> {

    let fp_reader = file_reader(input).map(fastq::Reader::new)?;
    if let Some(file) = input {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("skip first {} records", skip);
    info!("get {} records", take);

    let mut fp_writer =
        file_writer(output, compression_level, stdout_type).map(fastq::Writer::new)?;
    for rec in fp_reader
        .records()
        .skip(skip)
        .take(take)
        .map_while(Result::ok)
    {
        fp_writer.write_record(&rec)?;
    }
    fp_writer.flush()?;

    Ok(())
}
