use crate::utils::*;
use bio::io::fastq;
use anyhow::{Error, Ok};
use log::*;
use std::time::Instant;
use std::io::BufRead;

pub fn remove_read(
    file: &Option<&str>,
    out: &Option<&str>,
    name: &str,
) -> Result<(),Error> {
    info!("reading reads from file: {}", file.unwrap());
    info!("reading reads id form file: {}", name);
    let start = Instant::now();

    let mut ids = vec![];
    let list = file_reader(&Some(name))?;
    for i in list.lines().flatten(){
        ids.push(i);
    }
    
    if ids.len() == 0 {
        error!("reads id list is empty");
        std::process::exit(1);
    }

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut fq_writer = fastq::Writer::new(file_writer(out)?);
    for rec in fq_reader.records().flatten() {
        if !ids.contains(&rec.id().to_string()) {
            fq_writer.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;    
        }    
    }

    info!("time elapsed is: {:?}",start.elapsed());
    Ok(())
}