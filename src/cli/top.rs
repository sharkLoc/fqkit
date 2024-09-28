use crate::utils::*;
use anyhow::Result;
use bio::io::fastq;
use log::*;

pub fn top_n_records(
    input: Option<&String>,
    number: usize,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<()> {

    let fp = fastq::Reader::new(file_reader(input)?);
    if let Some(file) = input {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("get top {} records", number);

    let mut fo = fastq::Writer::new(file_writer(output, compression_level, stdout_type)?);
    for rec in fp.records().take(number).map_while(Result::ok) {
        fo.write_record(&rec)?;
    }
    fo.flush()?;

    Ok(())
}
