use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use std::time::Instant;

#[allow(clippy::too_many_arguments)]
pub fn flatten_fq(
    file: Option<&String>,
    out: Option<&String>,
    flag: u8,
    sep: char,
    gap: bool,
    len: bool,
    gc: bool,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {
    let start = Instant::now();

    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }
    info!("flag value is: {}", flag);

    if flag == 0 || flag > 15 {
        error!("error flag numer: {}, flag range [1..15]", flag);
        std::process::exit(1);
    }

    let fields = get_flag(flag);
    let mut out_writer = file_writer(out, compression_level, stdout_type)?;

    for rec in fq_reader.records().map_while(Result::ok) {
        let read = [rec.id().as_bytes(), rec.seq(), "+".as_bytes(), rec.qual()];
        let res = fields.iter().map(|idx| read[*idx]).collect::<Vec<&[u8]>>();

        let mut out = Vec::new();
        for x in res {
            out.push(std::str::from_utf8(x)?.to_string());
        }
        if gap {
            out.push(rec.seq().iter().filter(|x| *x == &b'N').count().to_string());
        }
        if len {
            out.push(rec.seq().len().to_string());
        }
        if gc {
            let gc_count = rec
                .seq()
                .iter()
                .filter(|x| *x == &b'G' || *x == &b'C')
                .count();
            let gc_ratio = format!("{:.2}", gc_count as f64 / rec.seq().len() as f64 * 100.0);
            out.push(gc_ratio);
        }
        out_writer.write_all(out.join(sep.to_string().as_str()).as_bytes())?;
        out_writer.write_all("\n".as_bytes())?;
    }
    out_writer.flush()?;

    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}

fn get_flag(num: u8) -> Vec<usize> {
    let flags = format!("{:b}", num).chars().rev().collect::<Vec<char>>();

    let mut fields = vec![];
    for (i, k) in flags.iter().enumerate() {
        if k == &'1' {
            fields.push(i);
        }
    }
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    // 1, 2, 4, 8
    #[test]
    fn flag3() {
        // 3 = 1 + 2, index: [0,1]
        assert_eq!(get_flag(3), vec![0, 1]);
    }

    #[test]
    fn flag7() {
        assert_eq!(get_flag(7), vec![0, 1, 2]);
    }
}
