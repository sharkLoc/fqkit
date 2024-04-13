use crate::utils::*;
use anyhow::{Error, Ok, Result};
use bio::io::{fasta, fastq};
use log::*;
use std::collections::HashMap;
use std::time::Instant;

pub fn cut_adapter(
    input: Option<&String>,
    seqfile: &String,
    left: bool,
    miss: usize,
    out: Option<&String>,
    compression_level: u32,
) -> Result<(), Error> {
    let start = Instant::now();
    if let Some(file) = input {
        info!("reading seq from file: {}", seqfile);
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let mut seqs = HashMap::new();
    let seqfile_reader = file_reader(Some(seqfile)).map(fasta::Reader::new)?;

    let mut iters = seqfile_reader.records();
    while let Some(each) = iters.next() {
        let rec = each?;
        if seqs.contains_key(&rec.id().to_owned()) {
            warn!("found duplicate sequence id: {}, keep first one", rec.id());
            continue;
        } else {
            seqs.entry(rec.id().to_owned())
                .or_insert(rec.seq().to_owned());
        }
    }
    if seqs.is_empty() {
        error!("file {} is empty", seqfile);
        std::process::exit(1);
    }

    let fq_reader = file_reader(input).map(fastq::Reader::new)?;
    let mut fq_writer = file_writer(out, compression_level).map(fastq::Writer::new)?;
    let mut flag = false;
    for rec in fq_reader.records().map_while(Result::ok) {
        let read_len = rec.seq().len();
        for (_, seq) in seqs.iter() {
            let pat = seq.as_slice();

            if read_len >= pat.len() {
                if left {
                    let hanming = pat
                        .iter()
                        .zip(rec.seq().iter())
                        .enumerate()
                        .filter(|(_, (y, z))| y != z)
                        .count();
                    if hanming <= miss {
                        fq_writer.write(
                            rec.id(),
                            rec.desc(),
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
                        fq_writer.write(
                            rec.id(),
                            rec.desc(),
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
            fq_writer.write_record(&rec)?;
        }
    }
    fq_writer.flush()?;

    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
