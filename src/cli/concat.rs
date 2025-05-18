use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{error, info};
use paraseq::fastq;
use std::io::BufRead;

pub fn concat_fqstq_lane(
    r1_list: &String,
    r2_list: &String,
    out_r1: &String,
    out_r2: &String,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut vec1 = vec![];
    let mut vec2 = vec![];
    let fp1 = file_reader(Some(r1_list))?;
    let fp2 = file_reader(Some(r2_list))?;

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

    let mut out_writer1 = file_writer(Some(out_r1), compression_level, stdout_type)?;
    let mut out_writer2 = file_writer(Some(out_r2), compression_level, stdout_type)?;
    let mut pe_read = 0;

    for pe in vec1.iter().zip(vec2.iter()) {
        info!("concat pe reads from file {} and {}", pe.0, pe.1);
        let mut fq_reader1 = file_reader(Some(pe.0)).map(fastq::Reader::new)?;
        let mut fq_reader2 = file_reader(Some(pe.1)).map(fastq::Reader::new)?;
        let mut rset1 = fastq::RecordSet::default();
        let mut rset2 = fastq::RecordSet::default();

        while rset1.fill(&mut fq_reader1)? && rset2.fill(&mut fq_reader2)? {
            for (rec1, rec2) in rset1
                .iter()
                .map_while(Result::ok)
                .zip(rset2.iter().map_while(Result::ok))
            {
                pe_read += 1;
                out_writer1.write_all(rec1.id())?;
                out_writer1.write_all(b"\n")?;
                out_writer1.write_all(rec1.seq())?;
                out_writer1.write_all(b"\n+\n")?;
                out_writer1.write_all(rec1.qual())?;
                out_writer1.write_all(b"\n")?;

                out_writer2.write_all(rec2.id())?;
                out_writer2.write_all(b"\n")?;
                out_writer2.write_all(rec2.seq())?;
                out_writer2.write_all(b"\n+\n")?;
                out_writer2.write_all(rec2.qual())?;
                out_writer2.write_all(b"\n")?;
            }
        }
        out_writer1.flush()?;
        out_writer2.flush()?;
    }

    info!(
        "total concat {} PE reads from {} different lane files",
        pe_read,
        vec1.len()
    );

    Ok(())
}
