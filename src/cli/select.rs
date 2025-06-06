use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::info;
use std::collections::HashMap;

pub fn select_pe_fastq(
    fq1: &String,
    fq2: &String,
    out_r1: &String,
    out_r2: &String,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut read1_id = HashMap::new();
    let mut read2_id = HashMap::new();
    let fq_reader1 = file_reader(Some(fq1)).map(fastq::Reader::new)?;
    let fq_reader2 = file_reader(Some(fq2)).map(fastq::Reader::new)?;

    for rec in fq_reader1.records().map_while(Result::ok) {
        let k = rec.id().to_owned();
        read1_id.entry(k).or_insert(());
    }
    for rec in fq_reader2.records().map_while(Result::ok) {
        let k = rec.id().to_owned();
        read2_id.entry(k).or_insert(());
    }
    info!("output selected read1 file: {}", out_r1);
    info!("output selected read2 file: {}", out_r2);
    let mut out_writer1 =
        file_writer(Some(out_r1), compression_level, stdout_type).map(fastq::Writer::new)?;
    let mut out_writer2 =
        file_writer(Some(out_r2), compression_level, stdout_type).map(fastq::Writer::new)?;
    let (mut pe_r1, mut pe_r2) = (0usize, 0usize);

    let fq_reader1 = file_reader(Some(fq1)).map(fastq::Reader::new)?;
    for rec in fq_reader1.records().map_while(Result::ok) {
        if read1_id.contains_key(rec.id()) && read2_id.contains_key(rec.id()) {
            pe_r1 += 1;
            out_writer1.write_record(&rec)?;
        }
    }
    out_writer1.flush()?;

    let fq_reader2 = file_reader(Some(fq2)).map(fastq::Reader::new)?;
    for rec in fq_reader2.records().map_while(Result::ok) {
        if read2_id.contains_key(rec.id()) && read1_id.contains_key(rec.id()) {
            pe_r2 += 1;
            out_writer2.write_record(&rec)?;
        }
    }
    out_writer2.flush()?;
    assert_eq!(pe_r1, pe_r2);

    info!("total selected pe reads: {}", pe_r1);

    Ok(())
}
