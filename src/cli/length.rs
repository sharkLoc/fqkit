use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use std::collections::HashMap;
use std::time::Instant;

pub fn fq_length(
    file: Option<&String>,
    rev: bool,
    out: Option<&String>,
    compression_level: u32,
) -> Result<(), Error> {
    let start = Instant::now();

    let mut reads_len = HashMap::new();
    let mut total = 0usize;
    let fp_reader = file_reader(file).map(fastq::Reader::new)?;
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let mut fo = file_writer(out, compression_level)?;
    for rec in fp_reader.records().map_while(Result::ok) {
        let rlen = rec.seq().len();
        *reads_len.entry(rlen).or_insert(0usize) += 1;
        total += 1;
    }

    let mut sort_len: Vec<(&usize, &usize)> = reads_len.iter().collect();
    if rev {
        sort_len.sort_by(|x, y| y.0.cmp(x.0));
    } else {
        sort_len.sort_by_key(|x| x.0);
    }

    fo.write_all("lenth\tcount\n".as_bytes())?;
    for (k, v) in sort_len.iter() {
        fo.write_all(format!("{}\t{}\n", k, v).as_bytes())?;
    }
    fo.flush()?;
    info!("total scan reads number: {}", total);
    info!("time elapsed is: {:?}", start.elapsed());

    Ok(())
}
