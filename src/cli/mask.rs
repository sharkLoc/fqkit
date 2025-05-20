use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{debug, info};
use paraseq::{fastq, fastx::Record};

pub fn mask_fastq(
    file: Option<&String>,
    phred: u8,
    qual_limit: u8,
    nt: char,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let (mut mask_base, mut mask_read) = (0, 0);
    let mut fq_reader = file_reader(file).map(fastq::Reader::new)?;

    info!("low quality value： {}", qual_limit);
    info!("mask low quality bases with: {}", nt);

    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;
    let mut rset = fastq::RecordSet::default();

    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            let score_min = rec.qual().iter().min().ok_or(FqkitError::EmptyQualRecord)? - phred;
            if score_min > qual_limit {
                fq_writer.write_all(rec.id())?;
                fq_writer.write_all(b"\n")?;
                fq_writer.write_all(rec.seq())?;
                fq_writer.write_all(b"\n+\n")?;
                fq_writer.write_all(rec.qual())?;
                fq_writer.write_all(b"\n")?;
            } else {
                debug!("mask read record id: {}", rec.id_str());
                mask_read += 1;

                let mut seq = Vec::with_capacity(rec.seq().len());
                for (s, q) in rec.seq().iter().zip(rec.qual().iter()) {
                    if q - phred <= qual_limit {
                        seq.push(nt as u8);
                        mask_base += 1;
                    } else {
                        seq.push(*s);
                    }
                }
                fq_writer.write_all(rec.id())?;
                fq_writer.write_all(b"\n")?;
                fq_writer.write_all(seq.as_slice())?;
                fq_writer.write_all(b"\n+\n")?;
                fq_writer.write_all(rec.qual())?;
                fq_writer.write_all(b"\n")?;
            }
        }
    }
    fq_writer.flush()?;

    info!("total mask {} bases from {} reads", mask_base, mask_read);
    Ok(())
}
