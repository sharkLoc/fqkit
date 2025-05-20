use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
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
            fq_writer.write_all(rec1.id())?;
            fq_writer.write_all(b"\n")?;
            fq_writer.write_all(rec1.seq())?;
            fq_writer.write_all(b"\n+\n")?;
            fq_writer.write_all(rec1.qual())?;
            fq_writer.write_all(b"\n")?;

            fq_writer.write_all(rec2.id())?;
            fq_writer.write_all(b"\n")?;
            fq_writer.write_all(rec2.seq())?;
            fq_writer.write_all(b"\n+\n")?;
            fq_writer.write_all(rec2.qual())?;
            fq_writer.write_all(b"\n")?;
        }
    }

    fq_writer.flush()?;

    info!("total PE reads number: {}", num);
    Ok(())
}
