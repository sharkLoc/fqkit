use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{error, info};
use paraseq::{fastq, fastx::Record};
use rand::{Rng, prelude::*};
use rand_pcg::Pcg64;

// reduce much memory but cost more time
fn select_fastq(
    file: Option<&String>,
    n: usize,
    seed: u64,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut rng: rand_pcg::Lcg128Xsl64 = Pcg64::seed_from_u64(seed);
    let mut get: Vec<usize> = Vec::with_capacity(n);

    let mut fq_reader = fastq::Reader::new(file_reader(file)?);
    info!("rand seed: {}", seed);
    info!("subseq number: {}", n);
    info!("reduce much memory but cost more time");
    let mut rset = fastq::RecordSet::default();
    let mut order: usize = 0;
    while rset.fill(&mut fq_reader)? {
        for _ in rset.iter().map_while(Result::ok) {
            order += 1;
            if order < n {
                get.push(order);
            } else {
                let ret = rng.random_range(0..=order);
                if ret < n {
                    get[ret] = order;
                }
            }
        }
    }

    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;
    let mut fq_reader2 = fastq::Reader::new(file_reader(file)?);
    let mut rset2 = fastq::RecordSet::default();
    let mut order2: usize = 0;

    while rset2.fill(&mut fq_reader2)? {
        for rec in rset2.iter().map_while(Result::ok) {
            order2 += 1;
            if get.contains(&order2) {
                fq_writer.write_all(rec.id())?;
                fq_writer.write_all(b"\n")?;
                fq_writer.write_all(rec.sep())?;
                fq_writer.write_all(b"\n+\n")?;
                fq_writer.write_all(rec.qual())?;
                fq_writer.write_all(b"\n")?;
            }
        }
    }
    fq_writer.flush()?;

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
) -> Result<(), FqkitError> {
    info!("rand seed: {}", seed);
    info!("subseq num: {}", n);
    info!("fast mode but cost more memory");

    let mut rng = Pcg64::seed_from_u64(seed);
    let mut get = Vec::with_capacity(n);

    let mut fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut rset = fastq::RecordSet::default();
    let mut order: usize = 0;
    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            order += 1;
            if order < n {
                let rec_t= vec![rec.id_str().to_owned(),rec.seq_str().to_owned(), rec.qual_str().to_owned()];
                get.push(rec_t);
            } else {
                let ret = rng.random_range(0..=order);
                if ret < n {
                    get[ret] = vec![rec.id_str().to_owned(),rec.seq_str().to_owned(), rec.qual_str().to_owned()];
                }
            }
        }
    }

    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;
    for rec in get.iter() {
        fq_writer.write_all(rec[0].as_bytes())?;
        fq_writer.write_all(b"\n")?;
        fq_writer.write_all(rec[1].as_bytes())?;
        fq_writer.write_all(b"\n+\n")?;
        fq_writer.write_all(rec[2].as_bytes())?;
        fq_writer.write_all(b"\n")?;
    }
    fq_writer.flush()?;

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
) -> Result<(), FqkitError> {
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
