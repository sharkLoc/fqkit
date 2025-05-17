use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{error, info};
use paraseq::fastq;

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
) -> Result<(), FqkitError> {
    let mut fq_reader = file_reader(file).map(fastq::Reader::new)?;
    info!("flag value is: {}", flag);

    if flag == 0 || flag > 15 {
        error!("error flag numer: {}, flag range [1..15]", flag);
        std::process::exit(1);
    }

    let fields = get_flag(flag);
    let mut out_writer = file_writer(out, compression_level, stdout_type)?;
    let mut rset = fastq::RecordSet::default();
    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            let read = [rec.id(), rec.seq(), b"+", rec.qual()];
            let res: Vec<&[u8]> = fields.iter().map(|idx| read[*idx]).collect::<Vec<&[u8]>>();

            for x in res {
                out_writer.write_all(x)?;
                out_writer.write_all(sep.to_string().as_bytes())?;
            }
            if gap {
                let gc_count = rec.seq().iter().filter(|x| *x == &b'N').count();
                out_writer.write_all(format!("{}", gc_count).as_bytes())?;
                out_writer.write_all(sep.to_string().as_bytes())?;
            }
            if len {
                out_writer.write_all(format!("{}", rec.sep().len()).as_bytes())?;
                out_writer.write_all(sep.to_string().as_bytes())?;
            }
            if gc {
                let gc_count = rec
                    .seq()
                    .iter()
                    .filter(|x| *x == &b'G' || *x == &b'C')
                    .count();
                let gc_ratio = format!("{:.2}", gc_count as f64 / rec.seq().len() as f64 * 100.0);
                out_writer.write_all(gc_ratio.as_bytes())?;
            }
            out_writer.write_all(b"\n")?;
        }
    }
    out_writer.flush()?;

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
