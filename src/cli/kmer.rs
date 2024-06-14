use anyhow::Error;
use bio::io::fastq;
use log::*;
use std::{collections::HashMap, time::Instant};
use crate::utils::*;
use nthash::nthash;


pub fn kmer_count(
    input: Option<&String>,
    kmer_len: usize,
    header: bool,
    output: Option<&String>,
    compression_level: u32,
) -> Result<(),Error> {
    let start = Instant::now();
    let reader = file_reader(input).map(fastq::Reader::new)?;
    if let Some(file) = input {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let mut writer = file_writer(output, compression_level)?;
    let mut kmers = HashMap::new();
    let (mut sidx, mut eidx) = (0,kmer_len);

    for rec in reader.records().flatten() {
        let khash = nthash(rec.seq(), kmer_len);

        let end = rec.seq().len() - kmer_len + 1;
        while eidx <= end {
            let kseq = &rec.seq()[sidx..eidx];
            let khash_this = nthash(kseq, kmer_len)[0];
            for k in khash.iter() {
                if khash_this.eq(k) {
                    *kmers.entry(kseq.to_owned()).or_insert(0_u64) += 1;
                }
            }

            sidx += 1;
            eidx += 1;
        }
    }

    if header {
        writer.write_all("kmer\tcount\n".as_bytes())?;
    }
    for (k,v) in kmers {
        writer.write_all(k.as_slice())?;
        writer.write_all(format!("\t{}\n",v).as_bytes())?;
    }
    writer.flush()?;
    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}