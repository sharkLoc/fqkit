use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use rand::{prelude::*, Rng};
use rand_pcg::Pcg64;

// reduce much memory but cost more time
fn select_fastq(
    file: Option<&String>,
    n: usize,
    seed: u64,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {

    let mut rng = Pcg64::seed_from_u64(seed);
    let mut get: Vec<usize> = Vec::with_capacity(n);

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("rand seed: {}", seed);
    info!("subseq number: {}", n);
    info!("reduce much memory but cost more time");

    for (order, _) in fq_reader.records().map_while(Result::ok).enumerate() {
        if order < n {
            get.push(order);
        } else {
            let ret = rng.random_range(0..=order);
            if ret < n {
                get[ret] = order;
            }
        }
    }

    let fo = file_writer(out, compression_level, stdout_type)?;
    let mut w = fastq::Writer::new(fo);
    let fq_reader2 = fastq::Reader::new(file_reader(file)?);
    for (order, rec) in fq_reader2.records().map_while(Result::ok).enumerate() {
        if get.contains(&order) {
            w.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
        }
    }
    w.flush()?;

    Ok(())
}

// fast mode but cost more memory
fn select_fastq2(
    file: Option<&String>,
    n: usize,
    seed: u64,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("rand seed: {}", seed);
    info!("subseq num: {}", n);
    info!("fast mode but cost more memory");

    let mut rng = Pcg64::seed_from_u64(seed);
    let mut get: Vec<fastq::Record> = Vec::with_capacity(n);

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    for (order, rec) in fq_reader.records().map_while(Result::ok).enumerate() {
        if order < n {
            get.push(rec);
        } else {
            let ret = rng.random_range(0..=order);
            if ret < n {
                get[ret] = rec;
            }
        }
    }

    let fo = file_writer(out, compression_level, stdout_type)?;
    let mut w = fastq::Writer::new(fo);
    for rec in get {
        w.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
    }
    w.flush()?;

    Ok(())
}

pub fn subset_fastq(
    rdc: bool,
    file: Option<&String>,
    n: usize,
    seed: u64,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {
    if rdc {
        if file.is_none() {
            error!("opt -r used, fastq data can't from stdin.");
            std::process::exit(1);
        }
        select_fastq(file, n, seed, out, compression_level, stdout_type)?;
    } else {
        select_fastq2(file, n, seed, out, compression_level, stdout_type)?;
    }

    Ok(())
}
