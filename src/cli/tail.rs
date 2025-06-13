use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::info;
use paraseq::fastq;

pub fn tail_n_records(
    input: Option<&String>,
    number: usize,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut fq_reader = fastq::Reader::new(file_reader(input)?);
    info!("get tail {} records", number);

    let mut rset = fastq::RecordSet::default();
    let mut fq_writer = file_writer(output, compression_level, stdout_type)?;

    let mut total = 0usize;
    while rset.fill(&mut fq_reader)? {
        for _ in rset.iter().map_while(Result::ok) {
            total += 1;
        }
    }
    info!("fastq file total reads number: {}", total);

    let skip_n = if number >= total {
        total
    } else {
        total - number
    };
    let mut fq_reader2 = fastq::Reader::new(file_reader(input)?);
    let mut rset2 = fastq::RecordSet::default();
    let mut count = 0usize;
    while rset2.fill(&mut fq_reader2)? {
        for rec in rset2.iter().map_while(Result::ok) {
            if count >= skip_n {
                write_record(&mut fq_writer, rec.id(), rec.seq(), rec.qual())?;
            }
            count += 1;
        }
    }
    fq_writer.flush()?;

    Ok(())
}
