use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::{fasta, fastq};
use log::info;

pub fn fq2fa(
    file: Option<&String>,
    remove: bool,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
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
