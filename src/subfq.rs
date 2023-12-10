use crate::utils::*;
use bio::io::fastq;
use anyhow::{Error, Ok};
use log::*;
use rand::{prelude::*, Rng};
use rand_pcg::Pcg64;
use std::time::Instant;


// reduce much memory but cost more time
pub fn select_fastq(
    file: &Option<&str>, 
    n: usize, 
    seed: u64, 
    out: &Option<&str>,
    compression_level: u32,
    quiet: bool,
) -> Result<(),Error> {
    if !quiet {
        if let Some(file) = file {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
        info!("rand seed: {}",seed);
        info!("subseq number: {}", n);
        info!("reduce much memory but cost more time");
    }
    let start = Instant::now();

    let mut rng = Pcg64::seed_from_u64(seed);
    let mut get: Vec<usize> = Vec::with_capacity(n);

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    for (order, _) in fq_reader.records().flatten().enumerate() {
        if order < n {
            get.push(order);
        } else {
            let ret = rng.gen_range(0..=order);
            if ret < n {
                get[ret] = order;
            }
        }
    }

    let fo = file_writer(out, compression_level)?;
    let mut w = fastq::Writer::new(fo);
    let fq_reader2 = fastq::Reader::new(file_reader(file)?);
    for (order, rec) in fq_reader2.records().flatten().enumerate() {
        if get.contains(&order) {
            w.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
        }
    }
    w.flush()?;

    if !quiet{
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}

// fast mode but cost more memory
pub fn select_fastq2(
    file: &Option<&str>, 
    n: usize, 
    seed: u64, 
    out: &Option<&str>,
    compression_level: u32,
    quiet: bool,
) -> Result<(),Error> {
    if !quiet{
        if let Some(file) = file {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
        info!("rand seed: {}",seed);
        info!("subseq num: {}", n);
        info!("fast mode but cost more memory");
    }
    let start = Instant::now();

    let mut rng = Pcg64::seed_from_u64(seed);
    let mut get: Vec<fastq::Record> = Vec::with_capacity(n);

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    for (order, rec) in fq_reader.records().flatten().enumerate() {
        if order < n {
            get.push(rec);
        } else {
            let ret = rng.gen_range(0..=order);
            if ret < n {
                get[ret] = rec;
            }
        }
    }

    let fo = file_writer(out, compression_level)?;
    let mut w = fastq::Writer::new(fo);
    for rec in get {
        w.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
    }
    w.flush()?;
    
    if !quiet{
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}