use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::{fastq, fastq::Record};
use log::{error, info};

pub fn phred_score(
    file: Option<&String>,
    out: Option<&String>,
    to33: bool,
    to64: bool,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    if let Some(r) = file {
        info!("read file from: {}", r);
    } else {
        info!("read file from: stdin");
    }

    let mut n = 0;
    if to33 {
        n += 1;
    }
    if to64 {
        n += 1;
    }
    if n > 1 {
        error!("only one of the flags --to33 and --to64 is allowed");
        std::process::exit(1);
    }
    if n == 0 {
        error!("please specifiy one of the flags: --to33, --to64");
        std::process::exit(1);
    }

    let mut fq_writer = file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;
    for rec in fq_reader.records().map_while(Result::ok) {
        let mut qual = vec![];
        if to33 {
            for q in rec.qual() {
                qual.push(q - 31);
            }
        }
        if to64 {
            for q in rec.qual() {
                qual.push(q + 31);
            }
        }
        let record = Record::with_attrs(rec.id(), rec.desc(), rec.seq(), qual.as_slice());
        fq_writer.write_record(&record)?;
    }
    fq_writer.flush()?;

    Ok(())
}
