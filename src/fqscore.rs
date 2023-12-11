use crate::utils::*;
use bio::io::{fastq, fastq::Record};
use anyhow::{Result, Ok};
use log::*;
use std::time::Instant;


pub fn phred_score(
    file: &Option<&str>,
    out: &Option<&str>,
    to33: bool,
    to64: bool,
    compression_level: u32,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        if let Some(r) = file {
            info!("read file from: {}",r);
        } else {
            info!("read file from: stdin");
        }
    }
    
    if !quiet {
        let mut n = 0;
        if to33 { n += 1; }
        if to64 { n += 1; }
        if n > 1 {
            error!("only one of the flags --to33 and --to64 is allowed");
            std::process::exit(1);
        }
        if n == 0 {
            error!("please specifiy one of the flags: --to33, --to64");
            std::process::exit(1);
        }
    }

    let start = Instant::now();
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut fq_writer = file_writer(out, compression_level).map(fastq::Writer::new)?;

    for rec in fq_reader.records().flatten() {
        let mut qual = vec![];
        if to33 {
            for q in rec.qual() { qual.push(q-31);  }
        } 
        if to64 {
            for q in rec.qual() { qual.push(q+31);  }
        }
        let record = Record::with_attrs(rec.id(), rec.desc(), rec.seq(), qual.as_slice());
        fq_writer.write_record(&record)?;
    }
    fq_writer.flush()?;

    if !quiet {
        info!("time elapsed is: {:?}",start.elapsed());
    }

    Ok(())
}