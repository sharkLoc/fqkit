use crate::utils::*;
use bio::io::{fastq,fastq::Record};
use anyhow::{Ok, Result};
use log::*;
use std::time::Instant;


pub fn mask_fastq(
    file: Option<&String>,
    phred: u8,
    qual_limit: u8,
    nt: char,
    out: Option<&String>,
    compression_level: u32,
) -> Result<()> {
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("low quality value： {}", qual_limit);
    info!("mask low quality bases with: {}", nt);
    let start = Instant::now();

    let (mut mask_base, mut mask_read) = (0,0);
    let fp_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut fp_writer = file_writer(out, compression_level).map(fastq::Writer::new)?;
    
    for rec in fp_reader.records().flatten() {
        let score_min = rec.qual().iter().min().unwrap() - phred;
        if score_min > qual_limit {
            fp_writer.write_record(&rec)?;
        } else {
            trace!("mask read record id: {}", rec.id());
            mask_read += 1;

            let mut seq = String::new();
            for (s,q) in rec.seq().iter().zip(rec.qual().iter()) {
                if q - phred <= qual_limit {
                    seq.push(nt);
                    mask_base += 1;
                } else {
                    seq.push(*s as char);
                }
            }
            let record = Record::with_attrs(rec.id(), rec.desc(), seq.as_bytes(), rec.qual());
            fp_writer.write_record(&record)?;
        }
    }
    fp_writer.flush()?;

    
    info!("total mask {} bases from {} reads", mask_base, mask_read);
    info!("time elapsed is: {:?}",start.elapsed()); 
    Ok(())
}