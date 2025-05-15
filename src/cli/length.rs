use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::info;
use rayon::prelude::*;
use std::collections::HashMap;

pub fn fq_length(
    file: Option<&String>,
    rev: bool,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut reads_len = HashMap::new();
    let mut total = 0usize;
    let fp_reader = file_reader(file).map(fastq::Reader::new)?;

    let mut fo = file_writer(out, compression_level, stdout_type)?;
    for rec in fp_reader.records().map_while(Result::ok) {
        let rlen = rec.seq().len();
        *reads_len.entry(rlen).or_insert(0usize) += 1;
        total += 1;
    }

    let mut sort_len: Vec<(&usize, &usize)> = reads_len.iter().collect();
    if rev {
        sort_len.par_sort_by(|x, y| y.0.cmp(x.0));
    } else {
        sort_len.par_sort_by_key(|x| x.0);
    }

    fo.write_all("lenth\tcount\n".as_bytes())?;
    for (k, v) in sort_len.iter() {
        fo.write_all(format!("{}\t{}\n", k, v).as_bytes())?;
    }
    fo.flush()?;
    info!("total scan reads number: {}", total);

    Ok(())
}
