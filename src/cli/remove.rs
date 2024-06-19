use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use std::io::BufRead;
use std::time::Instant;

pub fn remove_read(
    file: Option<&String>,
    out: Option<&String>,
    name: &String,
    save: &String,
    rm: bool,
    compression_level: u32,
) -> Result<(), Error> {
    let start = Instant::now();

    let mut ids = vec![];
    let list = file_reader(Some(name))?;
    info!("reading reads id form file: {}", name);
    for i in list.lines().map_while(Result::ok) {
        ids.push(i);
    }
    if ids.is_empty() {
        error!("reads id list is empty");
        std::process::exit(1);
    }

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    if let Some(file) = file {
        info!("reading reads from file: {}", file);
    } else {
        info!("reading reads from stdin");
    }

    if !rm {
        info!("removed reads in file: {}", save);
    }

    let mut fq_writer = fastq::Writer::new(file_writer(out, compression_level)?);
    if rm {
        for rec in fq_reader.records().map_while(Result::ok) {
            if !ids.contains(&rec.id().to_string()) {
                fq_writer.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            }
        }
        fq_writer.flush()?;
    } else {
        let mut rm_writer = fastq::Writer::new(file_writer(Some(save), compression_level)?);
        for rec in fq_reader.records().map_while(Result::ok) {
            if !ids.contains(&rec.id().to_string()) {
                fq_writer.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            } else {
                rm_writer.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            }
        }
        fq_writer.flush()?;
        rm_writer.flush()?;
    }

    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
