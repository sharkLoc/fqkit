use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::info;
use paraseq::fastq;

pub fn range_fastq(
    input: Option<&String>,
    skip: usize,
    take: usize,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut fq_reader = file_reader(input).map(fastq::Reader::new)?;
    info!("skip first {} records", skip);
    info!("get {} records", take);

    let mut fq_writer = file_writer(output, compression_level, stdout_type)?;
    let mut rset = fastq::RecordSet::default();

    let mut skipped = 0;
    let mut taken = 0;

    'outer: while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            if skipped < skip {
                skipped += 1;
                continue;
            }
            if taken < take {
                write_record(&mut fq_writer, rec.id(), rec.seq(), rec.qual())?;
                taken += 1;
            }
            if taken >= take {
                break 'outer;
            }
        }
    }
    fq_writer.flush()?;
    Ok(())
}
