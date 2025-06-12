use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use super::misc::write_record;
use log::warn;
use paraseq::{fastq, fastx::Record};

pub fn trim_fq(
    file: Option<&String>,
    left: usize,
    right: usize,
    len: usize,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let length = right + left;
    let mut fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut rset = fastq::RecordSet::default();
    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;

    while rset.fill(&mut fq_reader)? {
        for (idx, rec) in rset.iter().map_while(Result::ok).enumerate() {
            let rlen = rec.seq().len();
            if left >= rlen || right >= rlen || length >= rlen {
                warn!(
                    "read: {} in order {} is short than {} , skip",
                    rec.id_str(),
                    idx + 1,
                    length
                );
                continue;
            }

            let end = rlen - right;
            let seq = &rec.seq()[left..end];
            let qual = &rec.qual()[left..end];
            if seq.len() < len {
                continue;
            }
            write_record(&mut fq_writer, rec.id(), seq, qual)?;
        }
    }
    // Flush the writer to ensure all data is written
    fq_writer.flush()?;

    Ok(())
}
