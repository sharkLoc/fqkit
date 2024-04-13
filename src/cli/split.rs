use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use std::time::Instant;

pub fn split_interleaved(
    file: Option<&String>,
    out_dir: &String,
    out_pre: &String,
    compression_level: u32,
) -> Result<(), Error> {
    let start = Instant::now();

    let pre1 = format!("{}/{}_r1.fq.gz", out_dir, out_pre);
    let pre2 = format!("{}/{}_r2.fq.gz", out_dir, out_pre);
    let mut fh1 = fastq::Writer::new(file_writer_append(&pre1, compression_level)?);
    let mut fh2 = fastq::Writer::new(file_writer_append(&pre2, compression_level)?);

    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("read1 output file: {}", pre1);
    info!("read2 output file: {}", pre2);

    let mut num = 0usize;
    let mut flag = true;
    let fq_reader = fastq::Reader::new(file_reader(file)?);
    for rec in fq_reader.records().flatten() {
        num += 1;
        if flag {
            fh1.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            flag = false;
        } else {
            fh2.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            flag = true;
        }
    }
    fh1.flush()?;
    fh2.flush()?;

    info!("total split PE reads number: {}", num);
    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
