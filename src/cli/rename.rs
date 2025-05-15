use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq::{self, Record};
use log::*;

#[allow(clippy::too_many_arguments)]
pub fn rename_fastq(
    input: Option<&String>,
    keep: bool,
    prefix: Option<String>,
    label: Option<&String>,
    before: bool,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let fp = fastq::Reader::new(file_reader(input)?);

    let mut fo = fastq::Writer::new(file_writer(output, compression_level, stdout_type)?);
    let mut n: usize = 0;

    let new_rec = |rec: Record, n: usize| -> Record {
        let newid = match &prefix {
            Some(pre) => {
                if before {
                    format!("{}{}{}", label.unwrap_or(&String::new()), pre, n)
                } else {
                    format!("{}{}{}", pre, label.unwrap_or(&String::new()), n)
                }
            }
            None => {
                if before {
                    format!("{}{}", label.unwrap_or(&String::new()), rec.id())
                } else {
                    format!("{}{}", rec.id(), label.unwrap_or(&String::new()))
                }
            }
        };

        if keep {
            Record::with_attrs(&newid, rec.desc(), rec.seq(), rec.qual())
        } else {
            Record::with_attrs(&newid, None, rec.seq(), rec.qual())
        }
    };

    for rec in fp.records().map_while(Result::ok) {
        n += 1;
        let record = new_rec(rec, n);
        fo.write_record(&record)?;
    }
    fo.flush()?;

    info!("total rename sequence number: {}", n);
    Ok(())
}
