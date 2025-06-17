use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use paraseq::fastq;
use super::misc::write_record;
use log::{error, info};
use std::io::BufRead;

pub fn remove_read(
    file: Option<&String>,
    out: Option<&String>,
    name: &String,
    save: &String,
    rm: bool,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut ids = vec![];
    let list = file_reader(Some(name))?;
    info!("reading reads id form file: {}", name);
    for i in list.lines().map_while(Result::ok) {
        ids.push(i.as_bytes().to_vec());
    }
    if ids.is_empty() {
        error!("reads id list is empty");
        std::process::exit(1);
    }

    let mut fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut rset = fastq::RecordSet::default();
    if !rm {
        info!("removed reads in file: {}", save);
    }

    let mut writer = file_writer(out, compression_level, stdout_type)?;
    if rm {
        while rset.fill(&mut fq_reader)? {
            for rec in rset.iter().map_while(Result::ok) {
                if !ids.contains(&rec.id().to_vec()) {
                    write_record(&mut writer, rec.id(), rec.seq(), rec.qual())?;
                }
            }
        }
        writer.flush()?;
    } else {
        let mut rm_writer = file_writer(Some(save), compression_level, stdout_type)?;
        while rset.fill(&mut fq_reader)? {
            for rec in rset.iter().map_while(Result::ok) {
                if !ids.contains(&rec.id().to_vec()) {
                    write_record(&mut writer, rec.id(), rec.seq(), rec.qual())?;
                } else {
                    write_record(&mut rm_writer, rec.id(), rec.seq(), rec.qual())?;
                }
            }
        }
        writer.flush()?;
        rm_writer.flush()?;
    }

    Ok(())
}
