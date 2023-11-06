use crate::utils::*;
use bio::io::fastq;
use anyhow::{Error, Ok};
use log::*;
use std::time::Instant;

pub fn split_interleaved(
    file: &Option<&str>,
    out_dir: &str,
    out_pre: &str,
    quiet: bool,
) -> Result<(),Error> {
    let start = Instant::now();

    let pre1 = format!("{}/{}_r1.fq.gz", out_dir, out_pre);
    let pre2 = format!("{}/{}_r2.fq.gz", out_dir, out_pre);
    let mut fh1 = fastq::Writer::new(file_writer_append(&pre1)?);
    let mut fh2 = fastq::Writer::new(file_writer_append(&pre2)?);

    if !quiet {
        info!("reading from file: {}", file.unwrap());
        info!("read1 output file: {}",pre1);
        info!("read2 output file: {}",pre2);
    }

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

    if !quiet {
        info!("total split PE reads number: {}", num);
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}
