use crate::utils::*;
use anyhow::{Result,Ok};
use bio::io::fastq;
use log::*;
use std::{
    time::Instant,
    collections::HashMap,
};


pub fn select_pe_fastq(
    fq1: &String,
    fq2: &String,
    out_r1: &String,
    out_r2: &String,
    compression_level: u32,
) -> Result<()> {
    info!("read forward reads from file: {}", fq1);
    info!("read reverse reads from file: {}", fq2);
    info!("output selected read1 file: {}", out_r1);
    info!("output selected read2 file: {}", out_r2);
    let start = Instant::now();

    let mut read1_id = HashMap::new();
    let mut read2_id = HashMap::new();
    let fq_reader1 = file_reader(Some(fq1)).map(fastq::Reader::new)?;
    let fq_reader2 = file_reader(Some(fq2)).map(fastq::Reader::new)?;

    for rec in fq_reader1.records().flatten() { 
        let k = rec.id().to_owned();
        read1_id.entry(k).or_insert(()); 
    }
    for rec in fq_reader2.records().flatten() { 
        let k = rec.id().to_owned();
        read2_id.entry(k).or_insert(()); 
    }
    

    let mut out_writer1 = file_writer(Some(out_r1), compression_level).map(fastq::Writer::new)?;
    let mut out_writer2 = file_writer(Some(out_r2), compression_level).map(fastq::Writer::new)?;
    let (mut pe_r1, mut pe_r2) = (0usize, 0usize);

    let fq_reader1 = file_reader(Some(fq1)).map(fastq::Reader::new)?;
    for rec in fq_reader1.records().flatten() {
        if read1_id.contains_key(rec.id()) && read2_id.contains_key(rec.id()) {
            pe_r1 += 1;
            out_writer1.write_record(&rec)?;
        }
    }
    out_writer1.flush()?;
    
    let fq_reader2 = file_reader(Some(fq2)).map(fastq::Reader::new)?;
    for rec in fq_reader2.records().flatten() {
        if read2_id.contains_key(rec.id()) && read1_id.contains_key(rec.id()) {
            pe_r2 += 1;
            out_writer2.write_record(&rec)?;
        }
    }
    out_writer2.flush()?;
    assert_eq!(pe_r1,pe_r2);

    info!("total selected pe reads: {}", pe_r1);
    info!("time elapsed is: {:?}",start.elapsed());

    Ok(())
}