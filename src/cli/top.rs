use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::info;
use paraseq::fastq;

pub fn top_n_records(
    input: Option<&String>,
    number: usize,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut fq_reader = fastq::Reader::new(file_reader(input)?);
    info!("get top {} records", number);

    let mut rset = fastq::RecordSet::default();
    let mut fq_writer = file_writer(output, compression_level, stdout_type)?;
    let mut count = 0usize;

    'outer: while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            if count >= number {
                break 'outer;
            }
            fq_writer.write_all(rec.id())?;
            fq_writer.write_all(b"\n")?;
            fq_writer.write_all(rec.seq())?;
            fq_writer.write_all(b"\n+\n")?;
            fq_writer.write_all(rec.qual())?;
            fq_writer.write_all(b"\n")?;
            count += 1;
        }
    }
    fq_writer.flush()?;

    Ok(())
}
