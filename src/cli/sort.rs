use super::misc::write_record;
use crate::{
    errors::FqkitError,
    utils::{file_reader, file_writer},
};
use paraseq::{fastq, fastx::Record};
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

    let mut fq_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut rset = fastq::RecordSet::default();
    let mut vec_reads = vec![];

    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            vec_reads.push(vec![
                rec.id_str().to_string(),
                rec.seq_str().to_string(),
                rec.qual_str().to_string(),
            ]);
        }
    }

    if vec_reads.is_empty() {
        error!("no records found in the input file");
        std::process::exit(1);
    }
    info!("all records has been readed into memory, start sort ...");
    if reverse {
        info!("output reversed result");
    }

    if sort_by_name {
        info!("sort reads by name");

        if reverse {
            vec_reads.par_sort_by(|a, b| b[0].cmp(&a[0]));
        } else {
            vec_reads.par_sort_by(|a, b| a[0].cmp(&b[0]));
        }
    } else if sort_by_seq {
        info!("sort reads by sequence");
        if reverse {
            vec_reads.par_sort_by(|a, b| b[1].cmp(&a[1]));
        } else {
            vec_reads.par_sort_by(|a, b| a[1].cmp(&b[1]));
        }
    } else if sort_by_length {
        info!("sort reads by length");
        if reverse {
            vec_reads.par_sort_by_key(|b| std::cmp::Reverse(b[1].len()));
        } else {
            vec_reads.par_sort_by_key(|a| a[1].len());
        }
    } else if sort_by_gc {
        info!("sort reads by GC content");
        if reverse {
            vec_reads.par_sort_by(|a, b| {
                let r1_gc = a[1].chars().filter(|&x| x == 'G' || x == 'C').count() as f64
                    / a[1].len() as f64;
                let r2_gc = b[1].chars().filter(|&x| x == 'G' || x == 'C').count() as f64
                    / b[1].len() as f64;
                r2_gc.partial_cmp(&r1_gc).unwrap()
            });
        } else {
            vec_reads.par_sort_by(|a, b| {
                let r1_gc = a[1].chars().filter(|&x| x == 'G' || x == 'C').count() as f64
                    / a[1].len() as f64;
                let r2_gc = b[1].chars().filter(|&x| x == 'G' || x == 'C').count() as f64
                    / b[1].len() as f64;
                r1_gc.partial_cmp(&r2_gc).unwrap()
            });
        }
    }

    info!("sort done, start to output ...");
    let mut writer = file_writer(out, compression_level, stdout_type)?;
    for rec in vec_reads {
        write_record(
            &mut writer,
            rec[0].as_bytes(),
            rec[1].as_bytes(),
            rec[2].as_bytes(),
        )?;
    }
    writer.flush()?;

    Ok(())
}
