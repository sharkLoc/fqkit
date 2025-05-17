use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::info;
use paraseq::fastq;

pub fn fq2fa(
    file: Option<&String>,
    remove: bool,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut num = 0usize;
    let mut fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut rset = fastq::RecordSet::default();

    let mut fa_writer = file_writer(out, compression_level, stdout_type)?;

    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            num += 1;
            if remove {
                let id = rec.id().splitn(2, |&e| e == b' ').next().unwrap();
                fa_writer.write_all(id)?;
            } else {
                fa_writer.write_all(rec.id())?;
            }
            fa_writer.write_all(b"\n")?;
            fa_writer.write_all(rec.seq())?;
            fa_writer.write_all(b"\n")?;
        }
    }
    fa_writer.flush()?;

    info!("total reads number: {}", num);
    Ok(())
}
