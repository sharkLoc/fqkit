use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::info;
use std::collections::HashMap;

pub fn reverse_comp_seq(
    input: Option<&String>,
    out: Option<&String>,
    rev: bool,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let fq_reader = file_reader(input).map(fastq::Reader::new)?;
    if let Some(file) = input {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let maps = HashMap::from([
        (b'A', b'T'),
        (b'T', b'A'),
        (b'G', b'C'),
        (b'C', b'G'),
        (b'N', b'N'),
    ]);
    let mut out_writer =
        file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;

    for rec in fq_reader.records().map_while(Result::ok) {
        let rev_seq = rec.seq().iter().copied().rev().collect::<Vec<u8>>();
        let rev_qual = rec.qual().iter().copied().rev().collect::<Vec<u8>>();

        let rc_seq = rev_seq
            .iter()
            .map(|x| maps.get(x).unwrap_or(&b'N'))
            .collect::<Vec<&u8>>();
        let rev_comp = rc_seq.iter().map(|x| **x).collect::<Vec<u8>>();

        if rev {
            out_writer.write(
                rec.id(),
                rec.desc(),
                rev_seq.as_slice(),
                rev_qual.as_slice(),
            )?;
        } else {
            out_writer.write(
                rec.id(),
                rec.desc(),
                rev_comp.as_slice(),
                rev_qual.as_slice(),
            )?;
        }
    }
    out_writer.flush()?;

    Ok(())
}
