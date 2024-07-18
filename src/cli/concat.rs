use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use std::{io::BufRead, time::Instant};

pub fn concat_fqstq_lane(
    r1_list: &String,
    r2_list: &String,
    out_r1: &String,
    out_r2: &String,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {
    let start = Instant::now();

    let mut vec1 = vec![];
    let mut vec2 = vec![];
    let fp1 = file_reader(Some(r1_list))?;
    let fp2 = file_reader(Some(r2_list))?;
    info!("read forward reads for lane list: {}", r1_list);
    info!("read forward reads for lane list: {}", r2_list);

    for r1 in fp1.lines().map_while(Result::ok) {
        vec1.push(r1);
    }
    for r2 in fp2.lines().map_while(Result::ok) {
        vec2.push(r2);
    }
    if vec1.len() != vec2.len() {
        error!(
            "the number of fastq files in {} and {} is not equal",
            r1_list, r2_list
        );
        std::process::exit(1);
    }

    info!("outout read1 in file: {}", out_r1);
    info!("outout read1 in file: {}", out_r2);

    let mut out_writer1 =
        file_writer(Some(out_r1), compression_level, stdout_type).map(fastq::Writer::new)?;
    let mut out_writer2 =
        file_writer(Some(out_r2), compression_level, stdout_type).map(fastq::Writer::new)?;
    let mut pe_read = 0;
    for pe in vec1.iter().zip(vec2.iter()) {
        info!("concat pe reads from file {} and {}", pe.0, pe.1);
        let fq_reader1 = file_reader(Some(pe.0)).map(fastq::Reader::new)?;
        let fq_reader2 = file_reader(Some(pe.1)).map(fastq::Reader::new)?;
        for (rec1, rec2) in fq_reader1
            .records()
            .map_while(Result::ok)
            .zip(fq_reader2.records().map_while(Result::ok))
        {
            out_writer1.write_record(&rec1)?;
            out_writer2.write_record(&rec2)?;
            pe_read += 1;
        }
        out_writer1.flush()?;
        out_writer2.flush()?;
    }

    info!(
        "total concat {} PE reads from {} different lane files",
        pe_read,
        vec1.len()
    );
    info!("time elapsed is: {:?}", start.elapsed());

    Ok(())
}
