use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use std::{path::PathBuf, time::Instant};

pub fn split_chunk(
    file: Option<&String>,
    num: usize,
    gzip: bool,
    bzip2: bool,
    xz: bool,
    out_pre: &str,
    out_dir: &str,
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

    let (mut flag, mut index) = (0usize, 0usize);
    let out = if gzip {
        //format!("{}/{}{}.fastq.gz", out_dir, out_pre, index)
        PathBuf::from(out_dir).join(format!("{}{}.fq.gz",out_pre,index))
    } else if bzip2 {
        //format!("{}/{}{}.fastq.bz2", out_dir, out_pre, index)
        PathBuf::from(out_dir).join(format!("{}{}.fq.bz2",out_pre,index))
    } else if xz {
        //format!("{}/{}{}.fastq.xz", out_dir, out_pre, index)
        PathBuf::from(out_dir).join(format!("{}{}.fq.xz",out_pre,index))
    } else {
        //format!("{}/{}{}.fastq", out_dir, out_pre, index)
        PathBuf::from(out_dir).join(format!("{}{}.fq",out_pre,index))
    };

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    let mut fh = vec![fastq::Writer::new(file_writer(
        Some(&out),
        compression_level,
    )?)];

    info!("start to write file: {:?}", out);
    for rec in fq_reader.records().flatten() {
        if flag < num {
            let fhthis = fh.get_mut(index).unwrap();
            fhthis.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            flag += 1;
        } else {
            index += 1;
            let out = if gzip {
                PathBuf::from(out_dir).join(format!("{}{}.fq.gz",out_pre,index))
            } else if bzip2 {
                PathBuf::from(out_dir).join(format!("{}{}.fq.bz2",out_pre,index))
            } else if xz {
                PathBuf::from(out_dir).join(format!("{}{}.fq.xz",out_pre,index))
            } else {
                PathBuf::from(out_dir).join(format!("{}{}.fq",out_pre,index))
            };
            fh.push(fastq::Writer::new(file_writer(
                Some(&out),
                compression_level,
            )?));
            let fhthis = fh.get_mut(index).unwrap();

            info!("start to write file: {:?}", out);
            fhthis.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            flag = 1; // already write one record in this loop, flag add one
        }
    }

    info!("total chunk number is: {}", index + 1);
    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
