use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::info;
use paraseq::fastq;
use std::collections::HashSet;

pub fn select_pe_fastq(
    fq1: &String,
    fq2: &String,
    out_r1: &String,
    out_r2: &String,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut fq_reader1 = file_reader(Some(fq1)).map(fastq::Reader::new)?;
    let mut fq_reader2 = file_reader(Some(fq2)).map(fastq::Reader::new)?;
    let mut rset1 = fastq::RecordSet::default();
    let mut rset2 = fastq::RecordSet::default();

    let mut read1_id = HashSet::new();
    let mut read2_id = HashSet::new();
    while rset1.fill(&mut fq_reader1)? {
        for rec in rset1.iter().map_while(Result::ok) {
            read1_id.insert(rec.id().to_owned());
        }
    }

    while rset2.fill(&mut fq_reader2)? {
        for rec in rset2.iter().map_while(Result::ok) {
            read2_id.insert(rec.id().to_owned());
        }
    }

    let mut out_writer1 = file_writer(Some(out_r1), compression_level, stdout_type)?;
    let mut out_writer2 = file_writer(Some(out_r2), compression_level, stdout_type)?;

    let (mut pe_r1, mut pe_r2) = (0usize, 0usize);
    let intersect: HashSet<_> = read1_id.intersection(&read2_id).collect();

    let mut fq_reader1 = file_reader(Some(fq1)).map(fastq::Reader::new)?;
    info!("output selected read1 file: {}", out_r1);
    while rset1.fill(&mut fq_reader1)? {
        for rec in rset1.iter().map_while(Result::ok) {
            if intersect.contains(&rec.id().to_owned()) {
                pe_r1 += 1;
                write_record(&mut out_writer1, rec.id(), rec.seq(), rec.qual())?;
            }
        }
    }
    out_writer1.flush()?;

    let mut fq_reader2 = file_reader(Some(fq2)).map(fastq::Reader::new)?;
    info!("output selected read2 file: {}", out_r2);
    while rset2.fill(&mut fq_reader2)? {
        for rec in rset2.iter().map_while(Result::ok) {
            if intersect.contains(&rec.id().to_owned()) {
                pe_r2 += 1;
                write_record(&mut out_writer2, rec.id(), rec.seq(), rec.qual())?;
            }
        }
    }

    out_writer2.flush()?;
    assert_eq!(pe_r1, pe_r2);

    info!("total selected pe reads: {}", pe_r1);

    Ok(())
}
