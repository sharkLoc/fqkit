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
    save: &str,
    quiet: bool,
) -> Result<(),Error> {
    if !quiet {
        if let Some(file) = file {
            info!("reading reads from file: {}", file);
        } else {
            info!("reading reads from stdin");
        }
        info!("reading reads id form file: {}", name);
        info!("removed reads in file: {}", save);
    }
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
    let mut rm_writer = fastq::Writer::new(file_writer(&Some(&save))?);
    for rec in fq_reader.records().flatten() {
        if !ids.contains(&rec.id().to_string()) {
            fq_writer.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
        } else {
            rm_writer.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
        }    
    }

    if !quiet {
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}