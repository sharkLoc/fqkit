use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use regex::Regex;
use crate::utils::*;
use std::time::Instant;


pub fn search_fq(
    fq: &Option<&str>,
    pat: &str,
    out: &Option<&str>,
    quiet: bool,
) -> Result<(), Error> {
    let start = Instant::now();
    if !quiet {
        if let Some(file) = fq {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
        info!("regex pattern is: {}",pat);
        if let Some(out) = out {
            info!("reads write to file: {}", out);
        } else {
            info!("reads write to stdout");
        }
    }

    let re = Regex::new(pat)?;
    let fq_reader = file_reader(fq)
        .map(fastq::Reader::new)?;
    let mut fo = file_writer(out)
        .map(fastq::Writer::new)?;
    
    let mut num = 0usize;
    for rec in fq_reader.records().flatten() {
        let fq_str = std::str::from_utf8(rec.seq()).unwrap();
        if  re.is_match(fq_str) {
            num += 1;
            fo.write_record(&rec)?;
        }
    }

    if !quiet {
        info!("total reads matched number: {}",num);
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}