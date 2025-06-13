use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer_append};
use log::{error, info};
use paraseq::fastq;
use std::path::PathBuf;

pub fn split_interleaved(
    file: Option<&String>,
    out_dir: &String,
    out_pre: &String,
    gzip: bool,
    bzip2: bool,
    xz: bool,
    compression_level: u32,
) -> Result<(), FqkitError> {
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
    let mut rset = fastq::RecordSet::default();
    let mut fq_reader = fastq::Reader::new(file_reader(file)?);

    let pre1 = if gzip {
        PathBuf::from(out_dir).join(format!("{}_r1.fq.gz", out_pre))
    } else if bzip2 {
        PathBuf::from(out_dir).join(format!("{}_r1.fq.bz2", out_pre))
    } else if xz {
        PathBuf::from(out_dir).join(format!("{}_r1.fq.xz", out_pre))
    } else {
        PathBuf::from(out_dir).join(format!("{}_r1.fq", out_pre))
    };

    let pre2 = if gzip {
        PathBuf::from(out_dir).join(format!("{}_r2.fq.gz", out_pre))
    } else if bzip2 {
        PathBuf::from(out_dir).join(format!("{}_r2.fq.bz2", out_pre))
    } else if xz {
        PathBuf::from(out_dir).join(format!("{}_r2.fq.xz", out_pre))
    } else {
        PathBuf::from(out_dir).join(format!("{}_r2.fq", out_pre))
    };
    let mut fh1 = file_writer_append(&pre1, compression_level)?;
    info!("read1 output file: {}", pre1.display());
    let mut fh2 = file_writer_append(&pre2, compression_level)?;
    info!("read2 output file: {}", pre2.display());

    let mut num = 0usize;
    let mut flag = true;

    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            num += 1;
            if flag {
                write_record(&mut fh1, rec.id(), rec.seq(), rec.qual())?;
                flag = false;
            } else {
                write_record(&mut fh2, rec.id(), rec.seq(), rec.qual())?;
                flag = true;
            }
        }
    }
    fh1.flush()?;
    fh2.flush()?;

    info!("total split PE reads number: {}", num);
    Ok(())
}
