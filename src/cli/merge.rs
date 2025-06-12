use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use super::misc::write_record;
use log::info;
use paraseq::fastq;

pub fn interleaved(
    file1: &String,
    file2: &String,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut num: usize = 0usize;
    let mut fq1_reader = fastq::Reader::new(file_reader(Some(file1))?);
    let mut fq2_reader = fastq::Reader::new(file_reader(Some(file2))?);
    let mut rset1 = fastq::RecordSet::default();
    let mut rset2 = fastq::RecordSet::default();

    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;
    while rset1.fill(&mut fq1_reader)? && rset2.fill(&mut fq2_reader)? {
        for (rec1, rec2) in rset1
            .iter()
            .map_while(Result::ok)
            .zip(rset2.iter().map_while(Result::ok))
        {
            num += 2;
            write_record(&mut fq_writer, rec1.id(), rec1.seq(), rec1.qual())?;
            write_record(&mut fq_writer, rec2.id(), rec2.seq(), rec2.qual())?;
        }
    }

    fq_writer.flush()?;

    info!("total PE reads number: {}", num);
    Ok(())
}
