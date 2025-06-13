use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::*;
use paraseq::fastq;
use std::path::PathBuf;

#[allow(clippy::too_many_arguments)]
pub fn split_chunk(
    file: Option<&String>,
    num: usize,
    gzip: bool,
    bzip2: bool,
    xz: bool,
    out_pre: &str,
    out_dir: &str,
    compression_level: u32,
    stdout_type: char,
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

    let (mut flag, mut index) = (0usize, 0usize);
    let out = if gzip {
        PathBuf::from(out_dir).join(format!("{}{}.fq.gz", out_pre, index))
    } else if bzip2 {
        PathBuf::from(out_dir).join(format!("{}{}.fq.bz2", out_pre, index))
    } else if xz {
        PathBuf::from(out_dir).join(format!("{}{}.fq.xz", out_pre, index))
    } else {
        PathBuf::from(out_dir).join(format!("{}{}.fq", out_pre, index))
    };

    let mut fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut rest = fastq::RecordSet::default();
    let mut fh = vec![file_writer(Some(&out), compression_level, stdout_type)?];

    info!("start to write file: {}", out.display());
    while rest.fill(&mut fq_reader)? {
        for rec in rest.iter().map_while(Result::ok) {
            if flag < num {
                let mut this_writer = fh.get_mut(index).unwrap();
                write_record(&mut this_writer, rec.id(), rec.seq(), rec.qual())?;
                flag += 1;
            } else {
                index += 1;
                let out = if gzip {
                    PathBuf::from(out_dir).join(format!("{}{}.fq.gz", out_pre, index))
                } else if bzip2 {
                    PathBuf::from(out_dir).join(format!("{}{}.fq.bz2", out_pre, index))
                } else if xz {
                    PathBuf::from(out_dir).join(format!("{}{}.fq.xz", out_pre, index))
                } else {
                    PathBuf::from(out_dir).join(format!("{}{}.fq", out_pre, index))
                };

                fh.push(file_writer(Some(&out), compression_level, stdout_type)?);
                let mut fhthis = fh.get_mut(index).unwrap();

                info!("start to write file: {:?}", out);
                write_record(&mut fhthis, rec.id(), rec.seq(), rec.qual())?;
                flag = 1; // already write one record in this loop, flag add one
            }
        }
    }

    info!("total chunk number is: {}", index + 1);
    Ok(())
}
