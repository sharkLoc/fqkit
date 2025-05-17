use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::error;
use paraseq::fastq;

pub fn phred_score(
    file: Option<&String>,
    out: Option<&String>,
    to33: bool,
    to64: bool,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut fq_reader = file_reader(file).map(fastq::Reader::new)?;

    let mut n = 0;
    if to33 {
        n += 1;
    }
    if to64 {
        n += 1;
    }
    if n > 1 {
        error!("only one of the flags --to33 and --to64 is allowed");
        std::process::exit(1);
    }
    if n == 0 {
        error!("please specifiy one of the flags: --to33, --to64");
        std::process::exit(1);
    }

    let mut rset = fastq::RecordSet::default();
    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;

    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            let mut qual = vec![];
            if to33 {
                for q in rec.qual() {
                    qual.push(q - 31);
                }
            }
            if to64 {
                for q in rec.qual() {
                    qual.push(q + 31);
                }
            }
            fq_writer.write_all(rec.id())?;
            fq_writer.write_all(b"\n")?;
            fq_writer.write_all(rec.seq())?;
            fq_writer.write_all(b"\n+\n")?;
            fq_writer.write_all(&qual)?;
            fq_writer.write_all(b"\n")?;
        }
    }
    fq_writer.flush()?;

    Ok(())
}
