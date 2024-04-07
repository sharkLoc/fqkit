use anyhow::{Ok, Error};
use bio::io::fastq;
use std::collections::HashMap;
use log::*;
use crate::utils::*;
use std::time::Instant;


pub fn fq_length(
    file: Option<&String>,
    out: Option<&String>,
    compression_level: u32,
) -> Result<(),Error> {
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    let start = Instant::now();
    let mut reads_len = HashMap::new();
    let mut total = 0usize; 
    let fp_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut fo = file_writer(out, compression_level)?;

    for rec in fp_reader.records().flatten() {
        let rlen = rec.seq().len();
        *reads_len.entry(rlen).or_insert(0usize) += 1;
        total += 1;
    }

    let mut sort_len: Vec<(&usize,&usize)> = reads_len.iter().collect();
    sort_len.sort_by_key(|x| x.0);
    
    //fo.write("lenth\tcount\n".as_bytes())?;
    fo.write_all("lenth\tcount\n".as_bytes())?;
    for (k,v) in sort_len.iter() {
        fo.write_all(format!("{}\t{}\n",k,v).as_bytes())?;
    }
    fo.flush()?; 
    info!("total scan reads number: {}", total);
    info!("time elapsed is: {:?}", start.elapsed()); 

    Ok(())
}