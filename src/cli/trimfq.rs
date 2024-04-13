use crate::utils::*;
use anyhow::Result;
use bio::io::fastq;
use log::*;
use std::time::Instant;

pub fn trim_fq(
    file: Option<&String>,
    left: usize,
    right: usize,
    len: usize,
    out: Option<&String>,
    compression_level: u32,
) -> Result<()> {
    let start = Instant::now();
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let length = right + left;
    let fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut fq_writer = fastq::Writer::new(file_writer(out, compression_level)?);

    for (idx, rec) in fq_reader.records().flatten().enumerate() {
        let rlen = rec.seq().len();
        if left >= rlen || right >= rlen || length >= rlen {
            warn!(
                "read: {} in order {} is short than {} , skip",
                rec.id(),
                idx + 1,
                length
            );
            continue;
        }

        let end = rlen - right;
        let seq = &rec.seq()[left..end];
        let qual = &rec.qual()[left..end];
        if seq.len() < len {
            continue;
        }
        fq_writer.write(rec.id(), rec.desc(), seq, qual)?;
    }
    fq_writer.flush()?;

    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
