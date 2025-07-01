use crate::{
    cli::misc::{reverse_complement, write_record},
    errors::FqkitError,
    utils::{file_reader, file_writer_append},
};
use log::{error, info};
use paraseq::fastq;
use std::{
    collections::HashMap,
    io::BufRead,
    path::{Path, PathBuf},
};

fn barcode_list(file: &String, rev_comp: bool) -> Result<HashMap<Vec<u8>, String>, FqkitError> {
    let mut maps = HashMap::new();
    let fp = file_reader(Some(file))?;

    for line in fp.lines().map_while(Result::ok) {
        let item = line.split('\t').collect::<Vec<&str>>(); // barcode => sample
        let bar = if rev_comp {
            reverse_complement(item[0].as_bytes())
        } else {
            item[0]
                .chars()
                .rev()
                .collect::<String>()
                .as_bytes()
                .to_vec()
        };
        maps.entry(bar).or_insert(item[1].to_string());
    }
    Ok(maps)
}

#[inline]
fn hamming_dis(bar: &[u8], seq: &[u8]) -> usize {
    bar.iter().zip(seq.iter()).filter(|(a, b)| a != b).count()
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
) -> Result<(), FqkitError> {
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
            let (fq1, fq2, bar) = if gzip {
                (
                    PathBuf::from(outdir).join(format!("{}_1.fq.gz", name)),
                    PathBuf::from(outdir).join(format!("{}_2.fq.gz", name)),
                    PathBuf::from(outdir).join(format!("{}_barcode.fq.gz", name)),
                )
            } else if bzip2 {
                (
                    PathBuf::from(outdir).join(format!("{}_1.fq.bz2", name)),
                    PathBuf::from(outdir).join(format!("{}_2.fq.bz2", name)),
                    PathBuf::from(outdir).join(format!("{}_barcode.fq.bz2", name)),
                )
            } else if xz {
                (
                    PathBuf::from(outdir).join(format!("{}_1.fq.xz", name)),
                    PathBuf::from(outdir).join(format!("{}_2.fq.xz", name)),
                    PathBuf::from(outdir).join(format!("{}_barcode.fq.xz", name)),
                )
            } else {
                (
                    PathBuf::from(outdir).join(format!("{}_1.fq", name)),
                    PathBuf::from(outdir).join(format!("{}_2.fq", name)),
                    PathBuf::from(outdir).join(format!("{}_barcode.fq", name)),
                )
            };

            let fh1 = file_writer_append(&fq1, compression_level)?;
            let fh2 = file_writer_append(&fq2, compression_level)?;
            let fhb = file_writer_append(&bar, compression_level)?;
            fq_hand.push((bar_seq.clone(), bar_seq.len(), fh1, fh2, fhb));
        }

        info!("reading from read1 file: {}", big_fq1);
        let mut fq1_reader = fastq::Reader::new(file_reader(Some(big_fq1))?);
        info!("reading from read2 file: {}", big_fq2);
        let mut fq2_reader = fastq::Reader::new(file_reader(Some(big_fq2))?);
        let mut rset1 = fastq::RecordSet::default();
        let mut rset2 = fastq::RecordSet::default();
        let bar_count = fq_hand.len();
        let (mut read_pair, mut get_pair) = (0u64, 0u64);
        info!("barcode position mode: {}", mode);

        if mode == 2 {
            while rset1.fill(&mut fq1_reader)? && rset2.fill(&mut fq2_reader)? {
                for (rec1, rec2) in rset1
                    .iter()
                    .map_while(Result::ok)
                    .zip(rset2.iter().map_while(Result::ok))
                {
                    let read_len = rec2.seq().len();
                    read_pair += 1;

                    for idx in 0..bar_count {
                        if let Some((bar_seq, bar_len, fh1, fh2, fhb)) = fq_hand.get_mut(idx) {
                            if rec2.seq().len() < *bar_len {
                                continue;
                            }
                            let pos_index = read_len - *bar_len;
                            let read_part1 = &rec2.seq()[0..pos_index];
                            let read_part2 = &rec2.seq()[pos_index..];
                            let qual_part1 = &rec2.qual()[0..pos_index];
                            let qual_part2 = &rec2.qual()[pos_index..];
                            if hamming_dis(bar_seq, read_part2) <= mismatch {
                                get_pair += 1;
                                write_record(fh1, rec1.id(), rec1.seq(), rec1.qual())?;
                                write_record(fh2, rec2.id(), read_part1, qual_part1)?;
                                write_record(fhb, rec2.id(), read_part2, qual_part2)?;
                                break;
                            }
                        }
                    }
                }
            }
        } else if mode == 1 {
            while rset1.fill(&mut fq1_reader)? && rset2.fill(&mut fq2_reader)? {
                for (rec1, rec2) in rset1
                    .iter()
                    .map_while(Result::ok)
                    .zip(rset2.iter().map_while(Result::ok))
                {
                    read_pair += 1;
                    for idx in 0..bar_count {
                        if let Some((bar_seq, bar_len, fh1, fh2, fhb)) = fq_hand.get_mut(idx) {
                            if rec2.seq().len() < *bar_len {
                                continue;
                            }
                            let pos_index = *bar_len;
                            let read_part1 = &rec2.seq()[0..pos_index];
                            let read_part2 = &rec2.seq()[pos_index..];
                            let qual_part1 = &rec2.qual()[0..pos_index];
                            let qual_part2 = &rec2.qual()[pos_index..];
                            if hamming_dis(bar_seq, read_part1) <= mismatch {
                                get_pair += 1;
                                write_record(fh1, rec1.id(), rec1.seq(), rec1.qual())?;
                                write_record(fh2, rec2.id(), read_part2, qual_part2)?;
                                write_record(fhb, rec2.id(), read_part1, qual_part1)?;
                                break;
                            }
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
