use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fasta;
use log::{error, warn};
use paraseq::fastq;
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
    let seqfile_reader = file_reader(Some(seqfile)).map(fasta::Reader::new)?;

    let mut seqs = HashMap::new();
    let iters = seqfile_reader.records();
    for rec in iters.map_while(Result::ok) {
        if seqs.contains_key(rec.id()) {
            warn!("found duplicate sequence id: {}, keep first one", rec.id());
            continue;
        } else {
            seqs.entry(rec.id().to_owned())
                .or_insert(rec.seq().to_owned());
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
                            fq_writer.write_all(rec.id())?;
                            fq_writer.write_all(b"\n")?;
                            fq_writer.write_all(&rec.seq()[pat.len()..])?;
                            fq_writer.write_all(b"\n+\n")?;
                            fq_writer.write_all(&rec.qual()[pat.len()..])?;
                            fq_writer.write_all(b"\n")?;
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
                            fq_writer.write_all(rec.id())?;
                            fq_writer.write_all(b"\n")?;
                            fq_writer.write_all(&rec.seq()[0..idx])?;
                            fq_writer.write_all(b"\n+\n")?;
                            fq_writer.write_all(&rec.qual()[0..idx])?;
                            fq_writer.write_all(b"\n")?;
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
                fq_writer.write_all(rec.id())?;
                fq_writer.write_all(b"\n")?;
                fq_writer.write_all(rec.seq())?;
                fq_writer.write_all(b"\n+\n")?;
                fq_writer.write_all(rec.qual())?;
                fq_writer.write_all(b"\n")?;
            }
        }
    }
    fq_writer.flush()?;

    Ok(())
}
