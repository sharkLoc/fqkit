use crate::utils::*;
use bio::io::fastq;
use anyhow::{Error, Ok};
use log::*;
use std::time::Instant;


pub fn fq2fa(
    file: &Option<&str>, 
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
    let mut fo = file_writer(out)?;

    for rec in fq_reader.records().flatten() {
        num += 1;
        let pre = rec.id();
        let seq = std::str::from_utf8(rec.seq()).expect("Invalid UTF-8 sequence");
        let fa = match rec.desc() {
            Some(desc) => {
                format!(">{} {}\n{}\n", pre, desc, seq)
            }
            None => {
                format!(">{}\n{}\n", pre, seq)
            }
        };
        write!(&mut fo, "{}", fa)?;
    }

    if !quiet {
        info!("total reads number: {}",num);
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}