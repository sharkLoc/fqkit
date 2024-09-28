use crate::error::FqkitError;
use crate::utils::*;
use anyhow::Result;
use bio::io::fastq;
use log::*;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

fn barcode_list(file: &String, rev_comp: bool) -> Result<HashMap<String, String>> {
    let mut maps = HashMap::new();
    let mut error_flag = "";
    let fp = file_reader(Some(file))?;
    info!("reading from barcode list file: {}", file);

    if rev_comp {
        for line in fp.lines().map_while(Result::ok) {
            let item = line.split('\t').collect::<Vec<&str>>(); // barcode => sample
            let bar: String = item[0]
                .chars()
                .rev()
                .map(|x| match x {
                    'A' => 'T',
                    'G' => 'C',
                    'T' => 'A',
                    'C' => 'G',
                    'a' => 'T',
                    'g' => 'C',
                    't' => 'A',
                    'c' => 'G',
                    _ => {
                        warn!("invalid barcode base in sample {}, skipped", item[1]);
                        error_flag = "err";
                        'X'
                    }
                })
                .collect();
            if error_flag != "err" {
                maps.entry(bar).or_insert(item[1].to_string());
            }
            error_flag = "";
        }
    } else {
        for line in fp.lines().map_while(Result::ok) {
            let item = line.split('\t').collect::<Vec<&str>>();
            let bar: String = item[0]
                .chars()
                .map(|x| match x {
                    'A' => 'A',
                    'T' => 'T',
                    'G' => 'G',
                    'C' => 'C',
                    'a' => 'A',
                    't' => 'T',
                    'g' => 'G',
                    'c' => 'C',
                    _ => {
                        warn!("invalid barcode base in sample {}, skipped", item[1]);
                        error_flag = "err";
                        'X'
                    }
                })
                .collect();
            if error_flag != "err" {
                maps.entry(bar).or_insert(item[1].to_string());
            }
            error_flag = "";
        }
    }
    Ok(maps)
}

#[inline]
fn hamming_dis(bar: &str, seq: &[u8]) -> usize {
    let mut dist = 0usize;
    let bar = bar.as_bytes();
    for i in 0..bar.len() {
        if bar[i] != seq[i] {
            dist += 1;
        }
    }
    dist
}

#[allow(clippy::too_many_arguments)]
pub fn split_fq(
    big_fq1: &String,
    big_fq2: &String,
    bar_file: &String,
    rev_comp: bool,
    mode: usize,
    mismatch: usize,
    outdir: &str,
    gzip: bool,
    bzip2: bool,
    xz: bool,
    compression_level: u32,
) -> Result<()> {
    

    if !Path::new(outdir).try_exists().unwrap() {
        error!("{}", FqkitError::InvalidOutputDir(outdir.to_string()));
        std::process::exit(1);
    }
    let mut n = 0;
    if gzip {
        n += 1;
    }
    if bzip2 {
        n += 1;
    }
    if xz {
        n += 1;
    }
    if n > 1 {
        error!("only one of the flags --gzip, --xz and --bzip2 is allowed");
        std::process::exit(1);
    }

    if let Ok(maps) = barcode_list(bar_file, rev_comp) {
        if maps.is_empty() {
            error!("{}", FqkitError::EmptyFile(bar_file.to_string()));
            std::process::exit(1);
        }

        let mut fq_hand = Vec::new();
        for (bar_seq, name) in maps {
            let fq1 = if gzip {
                PathBuf::from(outdir).join(format!("{}_1.fq.gz", name))
            } else if bzip2 {
                PathBuf::from(outdir).join(format!("{}_1.fq.bz2", name))
            } else if xz {
                PathBuf::from(outdir).join(format!("{}_1.fq.xz", name))
            } else {
                PathBuf::from(outdir).join(format!("{}_1.fq", name))
            };
            let fq2 = if gzip {
                PathBuf::from(outdir).join(format!("{}_2.fq.gz", name))
            } else if bzip2 {
                PathBuf::from(outdir).join(format!("{}_2.fq.bz2", name))
            } else if xz {
                PathBuf::from(outdir).join(format!("{}_2.fq.xz", name))
            } else {
                PathBuf::from(outdir).join(format!("{}_2.fq", name))
            };
            let bar = if gzip {
                PathBuf::from(outdir).join(format!("{}_barcode.fq.gz", name))
            } else if bzip2 {
                PathBuf::from(outdir).join(format!("{}_barcode.fq.bz2", name))
            } else if xz {
                PathBuf::from(outdir).join(format!("{}_barcode.fq.xz", name))
            } else {
                PathBuf::from(outdir).join(format!("{}_barcode.fq", name))
            };

            let fh1 = fastq::Writer::new(file_writer_append(&fq1, compression_level)?);
            let fh2 = fastq::Writer::new(file_writer_append(&fq2, compression_level)?);
            let fhb = fastq::Writer::new(file_writer_append(&bar, compression_level)?);
            let len = bar_seq.len();
            fq_hand.push((bar_seq, len, fh1, fh2, fhb));
        }

        let bar_count = fq_hand.len();
        let fq1_reader = fastq::Reader::new(file_reader(Some(big_fq1))?);
        let fq2_reader = fastq::Reader::new(file_reader(Some(big_fq2))?);
        let (mut read_pair, mut get_pair) = (0u64, 0u64);

        info!("reading from read1 file: {}", big_fq1);
        info!("reading from read2 file: {}", big_fq2);
        info!("barcode position mode: {}", mode);

        if mode == 2 {
            for (rec1, rec2) in fq1_reader
                .records()
                .flatten()
                .zip(fq2_reader.records().flatten())
            {
                let read_len = rec2.seq().len();
                read_pair += 1;
                for idx in 0..bar_count {
                    if let Some((bar_seq, bar_len, fh1, fh2, fhb)) = fq_hand.get_mut(idx) {
                        let pos_index = read_len - *bar_len;
                        let read_part1 = &rec2.seq()[0..pos_index];
                        let read_part2 = &rec2.seq()[pos_index..];
                        let qual_part1 = &rec2.qual()[0..pos_index];
                        let qual_part2 = &rec2.qual()[pos_index..];
                        if hamming_dis(bar_seq, read_part2) <= mismatch {
                            get_pair += 1;
                            fh1.write(rec1.id(), rec1.desc(), rec1.seq(), rec1.qual())?;
                            fh2.write(rec2.id(), rec2.desc(), read_part1, qual_part1)?;
                            fhb.write(rec2.id(), rec2.desc(), read_part2, qual_part2)?;
                            break;
                        }
                    }
                }
            }
        } else if mode == 1 {
            for (rec1, rec2) in fq1_reader
                .records()
                .flatten()
                .zip(fq2_reader.records().flatten())
            {
                read_pair += 1;
                for idx in 0..bar_count {
                    if let Some((bar_seq, bar_len, fh1, fh2, fhb)) = fq_hand.get_mut(idx) {
                        let pos_index = *bar_len;
                        let read_part1 = &rec2.seq()[0..pos_index];
                        let read_part2 = &rec2.seq()[pos_index..];
                        let qual_part1 = &rec2.qual()[0..pos_index];
                        let qual_part2 = &rec2.qual()[pos_index..];
                        if hamming_dis(bar_seq, read_part1) <= mismatch {
                            get_pair += 1;
                            fh1.write(rec1.id(), rec1.desc(), rec1.seq(), rec1.qual())?;
                            fh2.write(rec2.id(), rec2.desc(), read_part2, qual_part2)?;
                            fhb.write(rec2.id(), rec2.desc(), read_part1, qual_part1)?;
                            break;
                        }
                    }
                }
            }
        } else {
            error!("invalid mode arg, must be 1 or 2 !");
            std::process::exit(1);
        }
        info!(
            "data split rate: {:.4}%",
            get_pair as f64 / read_pair as f64 * 100.0
        );
    }

    Ok(())
}
