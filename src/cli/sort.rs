use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::{error, info};
use rayon::prelude::*;

#[allow(clippy::too_many_arguments)]
pub fn sort_fastq(
    file: Option<&String>,
    sort_by_name: bool,
    sort_by_seq: bool,
    sort_by_gc: bool,
    sort_by_length: bool,
    reverse: bool,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut n = 0;
    if sort_by_gc {
        n += 1;
    }
    if sort_by_length {
        n += 1;
    }
    if sort_by_name {
        n += 1;
    }
    if sort_by_seq {
        n += 1;
    }
    if n > 1 {
        error!(
            "only one of the flags -l (--sort-by-length), -n (--sort-by-name), -g (--sort-by-gc) and -s (--sort-by-seq) is allowed"
        );
        std::process::exit(1);
    }
    if n == 0 {
        error!("please specifiy one of the flags: -l, -n, -g, -s");
        std::process::exit(1);
    }

    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    if reverse {
        info!("output reversed result");
    }

    let mut vec_reads = vec![];
    for rec in fq_reader.records().map_while(Result::ok) {
        vec_reads.push(rec);
    }
    info!("all records has been readed into memory, start sort ...");

    if sort_by_name {
        info!("sort read by name");
        if reverse {
            vec_reads.par_sort_by(|a, b| {
                let read_name1 = if let Some(des) = a.desc() {
                    format!("{} {}", a.id(), des)
                } else {
                    a.id().to_string()
                };
                let read_name2 = if let Some(des) = b.desc() {
                    format!("{} {}", b.id(), des)
                } else {
                    b.id().to_string()
                };
                read_name2.cmp(&read_name1)
            });
        } else {
            vec_reads.par_sort_by(|a, b| {
                let read_name1 = if let Some(des) = a.desc() {
                    format!("{} {}", a.id(), des)
                } else {
                    a.id().to_string()
                };
                let read_name2 = if let Some(des) = b.desc() {
                    format!("{} {}", b.id(), des)
                } else {
                    b.id().to_string()
                };
                read_name1.cmp(&read_name2)
            });
        }
    } else if sort_by_seq {
        info!("sort read by sequence");
        if reverse {
            vec_reads.par_sort_by(|a, b| b.seq().cmp(a.seq()));
        } else {
            vec_reads.par_sort_by(|a, b| a.seq().cmp(b.seq()));
        }
    } else if sort_by_length {
        info!("sort read by length");
        if reverse {
            //vec_reads.sort_by(|a, b| b.seq().len().cmp(&a.seq().len()));
            vec_reads.par_sort_by_key(|b| std::cmp::Reverse(b.seq().len()))
        } else {
            //vec_reads.sort_by(|a, b| a.seq().len().cmp(&b.seq().len()));
            vec_reads.par_sort_by_key(|a| a.seq().len())
        }
    } else if sort_by_gc {
        info!("sort read by gc content");
        if reverse {
            vec_reads.par_sort_by(|a, b| {
                let r1_gc = a
                    .seq()
                    .iter()
                    .filter(|x| x == &&b'G' || x == &&b'C')
                    .count() as f64
                    / a.seq().len() as f64;
                let r2_gc = b
                    .seq()
                    .iter()
                    .filter(|x| x == &&b'G' || x == &&b'C')
                    .count() as f64
                    / b.seq().len() as f64;
                r2_gc.partial_cmp(&r1_gc).unwrap()
            });
        } else {
            vec_reads.par_sort_by(|a, b| {
                let r1_gc = a
                    .seq()
                    .iter()
                    .filter(|x| x == &&b'G' || x == &&b'C')
                    .count() as f64
                    / a.seq().len() as f64;
                let r2_gc = b
                    .seq()
                    .iter()
                    .filter(|x| x == &&b'G' || x == &&b'C')
                    .count() as f64
                    / b.seq().len() as f64;
                r1_gc.partial_cmp(&r2_gc).unwrap()
            });
        }
    }

    info!("sort done, start to output ...");
    let mut fq_writer = file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;
    for rec in vec_reads {
        fq_writer.write_record(&rec)?;
    }
    fq_writer.flush()?;

    Ok(())
}
