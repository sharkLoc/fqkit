use crate::utils::*;
use bio::io::fastq;
use bio::io::fasta;
use anyhow::{Error, Ok};
use log::*;
use std::time::Instant;


pub fn fq2fa(
    file: &Option<&str>,
    remove: bool,
    out: &Option<&str>,
    quiet: bool,
) -> Result<(), Error> {
    if !quiet {
        if let Some(file) = file {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
    }   
    let start = Instant::now();
    let mut num = 0usize;
    
    let fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut fo = fasta::Writer::new(file_writer(out)?);
   
    if remove {
        for rec in fq_reader.records().flatten() {
            num += 1;
            fo.write(rec.id(), None, rec.seq())?;
        }
    } else {
        for rec in fq_reader.records().flatten() {
            num += 1;
            fo.write(rec.id(), rec.desc(), rec.seq())?;
        }
    }
    

    if !quiet {
        info!("total reads number: {}",num);
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}