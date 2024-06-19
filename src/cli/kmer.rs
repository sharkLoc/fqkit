use crate::utils::*;
use anyhow::Error;
use bio::io::fastq;
use log::*;
use nthash::nthash;
use std::{collections::HashMap, time::Instant};

pub fn kmer_count(
    input: Option<&String>,
    kmer_len: usize,
    header: bool,
    output: Option<&String>,
    compression_level: u32,
) -> Result<(), Error> {
    let start = Instant::now();
    let reader = file_reader(input).map(fastq::Reader::new)?;
    if let Some(file) = input {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let mut writer = file_writer(output, compression_level)?;
    let mut kmers = HashMap::new();

    for rec in reader.records().map_while(Result::ok) {
        let (mut sidx, mut eidx) = (0, kmer_len);
        let khash = nthash(rec.seq(), kmer_len);
        let len = rec.seq().len();

        while eidx <= len {
            let kseq = &rec.seq()[sidx..eidx];
            let khash_this = nthash(kseq, kmer_len)[0];
            if khash.contains(&khash_this) {
                *kmers.entry(kseq.to_owned()).or_insert(0_u64) += 1;
            }
            sidx += 1;
            eidx += 1;
        }
    }

    if header {
        writer.write_all("kmer\tcount\n".as_bytes())?;
    }
    for (k, v) in kmers {
        writer.write_all(k.as_slice())?;
        writer.write_all(format!("\t{}\n", v).as_bytes())?;
    }
    writer.flush()?;
    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
