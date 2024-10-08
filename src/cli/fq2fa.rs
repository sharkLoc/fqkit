use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fasta;
use bio::io::fastq;
use log::*;

pub fn fq2fa(
    file: Option<&String>,
    remove: bool,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {

    let mut num = 0usize;
    let fq_reader = fastq::Reader::new(file_reader(file)?);
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let mut fo = fasta::Writer::new(file_writer(out, compression_level, stdout_type)?);
    if remove {
        for rec in fq_reader.records().map_while(Result::ok) {
            num += 1;
            fo.write(rec.id(), None, rec.seq())?;
        }
    } else {
        for rec in fq_reader.records().flatten() {
            num += 1;
            fo.write(rec.id(), rec.desc(), rec.seq())?;
        }
    }
    fo.flush()?;

    info!("total reads number: {}", num);
    Ok(())
}
