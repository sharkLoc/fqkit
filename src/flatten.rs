use crate::utils::*;
use bio::io::fastq;
use anyhow::{Error, Ok};
use log::*;
use std::time::Instant;

pub fn flatten_fq(
    file: &Option<&str>, 
    out: &Option<&str>,
    flag: u8,
    sep: char,
    compression_level: u32,
    quiet: bool,
) -> Result<(),Error> {
    let start = Instant::now();
    if !quiet {
        if let Some(file) = file {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
        info!("flag value is: {}", flag);
    }

    if flag == 0 || flag>15 {
        error!("error flag numer: {}, flag range [1..15]",flag);
        std::process::exit(1);
    }
    
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut out_writer = file_writer(out, compression_level)?;
    let flags = format!("{:b}",flag).chars().rev().collect::<Vec<char>>();
    let mut fields = vec![];
    for (i,k) in flags.iter().enumerate() {
        if k == &'1' {
            fields.push(i);
        }
    }

    for rec in fq_reader.records().flatten() {
        let read = vec![rec.id().as_bytes(), rec.seq(), "+".as_bytes(), rec.qual()];
        let res = fields.iter().map(|idx| read[*idx]).collect::<Vec<&[u8]>>();
        let mut out = Vec::new();
        for x in res {
            out.push(std::str::from_utf8(x)?.to_string());
        }
        out_writer.write(out.join(sep.to_string().as_str()).as_bytes())?;
        out_writer.write("\n".as_bytes())?;
    }
    out_writer.flush()?;
    
    if !quiet{
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}