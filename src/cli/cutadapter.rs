use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{error, warn};
use paraseq::{fasta, fastq};
use std::collections::HashMap;

pub fn cut_adapter(
    input: Option<&String>,
    seqfile: &String,
    left: bool,
    miss: usize,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut seqfile_reader = file_reader(Some(seqfile)).map(fasta::Reader::new)?;
    let mut faset = fasta::RecordSet::default();
    let mut seqs = HashMap::new();

    while faset.fill(&mut seqfile_reader)? {
        for rec in faset.iter().map_while(Result::ok) {
            if seqs.contains_key(rec.id()) {
                warn!(
                    "found duplicate sequence id: {}, keep first one",
                    std::str::from_utf8(rec.id())?
                );
                continue;
            } else {
                seqs.entry(rec.id().to_owned())
                    .or_insert(rec.seq().to_vec());
            }
        }
    }

    if seqs.is_empty() {
        error!("{}", FqkitError::EmptyFile(seqfile.to_string()));
        std::process::exit(1);
    }

    let mut fq_reader = file_reader(input).map(fastq::Reader::new)?;
    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;
    let mut flag = false;
    let mut rset = fastq::RecordSet::default();

    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            let read_len = rec.seq().len();

            for (_, seq) in seqs.iter() {
                let pat = seq;
                if read_len >= pat.len() {
                    if left {
                        let hanming = pat
                            .iter()
                            .zip(rec.seq().iter())
                            .enumerate()
                            .filter(|(_, (y, z))| y != z)
                            .count();
                        if hanming <= miss {
                            write_record(
                                &mut fq_writer,
                                rec.id(),
                                &rec.seq()[pat.len()..],
                                &rec.qual()[pat.len()..],
                            )?;
                            flag = true;
                            break;
                        }
                    } else {
                        let idx = read_len - pat.len();
                        let hanming = pat
                            .iter()
                            .zip(rec.seq()[idx..].iter())
                            .enumerate()
                            .filter(|(_, (y, z))| y != z)
                            .count();
                        if hanming <= miss {
                            write_record(
                                &mut fq_writer,
                                rec.id(),
                                &rec.seq()[0..idx],
                                &rec.qual()[0..idx],
                            )?;
                            flag = true;
                            break;
                        }
                    }
                }
            }
            if flag {
                flag = false;
                continue;
            } else {
                write_record(&mut fq_writer, rec.id(), rec.seq(), rec.qual())?;
            }
        }
    }
    fq_writer.flush()?;

    Ok(())
}
