use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use regex::Regex;
use crate::utils::*;
use std::time::Instant;


pub fn search_fq(
    fq: &str,
    pat: &str,
    out: &Option<&str>
) -> Result<(), Error> {
    let start = Instant::now();
    info!("reading frim file: {}",fq);
    info!("regex pattern is: {}",pat);
    if let Some(out) = out {
        info!("reads write to file: {}", out);
    } else {
        info!("reads write to stdout");
    }

    let re = Regex::new(pat)?;
    let fq_reader = file_reader(&Some(fq))
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

    info!("total reads matched number: {}",num);
    info!("time elapsed is: {:?}",start.elapsed());
    Ok(())
}