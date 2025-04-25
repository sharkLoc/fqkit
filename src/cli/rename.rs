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
    if let Some(file) = input {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let mut fo = fastq::Writer::new(file_writer(output, compression_level, stdout_type)?);
    let mut n: usize = 0;

    if let Some(pre) = prefix {
        for rec in fp.records().map_while(Result::ok) {
            n += 1;
            /*let newid = match label {
                Some(x) => {
                    if before {
                        format!("{}{}{}",x,pre,n)
                    } else {
                        format!("{}{}{}",pre,x,n)
                    }
                },
                None => { format!("{}{}",pre,n) }
            };*/
            let newid = if before {
                format!("{}{}{}", label.unwrap_or(&String::new()), pre, n)
            } else {
                format!("{}{}{}", pre, label.unwrap_or(&String::new()), n)
            };

            let record = if keep {
                Record::with_attrs(&newid, rec.desc(), rec.seq(), rec.qual())
            } else {
                Record::with_attrs(&newid, None, rec.seq(), rec.qual())
            };
            fo.write_record(&record)?;
        }
        fo.flush()?;
    } else {
        for rec in fp.records().map_while(Result::ok) {
            n += 1;
            /*let newid = if let Some(x) = label {
                if before {
                    format!("{}{}",x,rec.id())
                } else {
                    format!("{}{}",rec.id(),x)
                }
            } else {
                format!("{}",rec.id())
            };*/
            let newid = if before {
                format!("{}{}", label.unwrap_or(&String::new()), rec.id())
            } else {
                format!("{}{}", rec.id(), label.unwrap_or(&String::new()))
            };
            let record = if keep {
                Record::with_attrs(&newid, rec.desc(), rec.seq(), rec.qual())
            } else {
                Record::with_attrs(&newid, None, rec.seq(), rec.qual())
            };
            fo.write_record(&record)?;
        }
        fo.flush()?;
    }

    info!("total rename sequence number: {}", n);
    Ok(())
}
