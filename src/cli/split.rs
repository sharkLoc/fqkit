use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use std::{path::PathBuf, time::Instant};

pub fn split_interleaved(
    file: Option<&String>,
    out_dir: &String,
    out_pre: &String,
    gzip: bool,
    bzip2: bool,
    xz: bool,
    compression_level: u32,
) -> Result<(), Error> {
    let start = Instant::now();

    let mut n = 0;
    if gzip {
        n += 1;
    }
    if bzip2 {
        n += 1;
    }
    if xz {
        n += 1;
    }
    if n > 1 {
        error!("only one of the flags --gzip , --xz and --bzip2 is allowed");
        std::process::exit(1);
    }
    let fq_reader = fastq::Reader::new(file_reader(file)?);
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    //let pre1 = format!("{}/{}_r1.fq.gz", out_dir, out_pre);
    let pre1 = if gzip {
        PathBuf::from(out_dir).join(format!("{}_r1.fq.gz", out_pre))
    } else if bzip2 {
        PathBuf::from(out_dir).join(format!("{}_r1.fq.bz2", out_pre))
    } else if xz {
        PathBuf::from(out_dir).join(format!("{}_r1.fq.xz", out_pre))
    } else {
        PathBuf::from(out_dir).join(format!("{}_r1.fq", out_pre))
    };
    //let pre2 = format!("{}/{}_r2.fq.gz", out_dir, out_pre);
    let pre2 = if gzip {
        PathBuf::from(out_dir).join(format!("{}_r2.fq.gz", out_pre))
    } else if bzip2 {
        PathBuf::from(out_dir).join(format!("{}_r2.fq.bz2", out_pre))
    } else if xz {
        PathBuf::from(out_dir).join(format!("{}_r2.fq.xz", out_pre))
    } else {
        PathBuf::from(out_dir).join(format!("{}_r2.fq", out_pre))
    };
    let mut fh1 = fastq::Writer::new(file_writer_append(&pre1, compression_level)?);
    info!("read1 output file: {:?}", pre1);
    let mut fh2 = fastq::Writer::new(file_writer_append(&pre2, compression_level)?);
    info!("read2 output file: {:?}", pre2);

    let mut num = 0usize;
    let mut flag = true;
    for rec in fq_reader.records().map_while(Result::ok) {
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
