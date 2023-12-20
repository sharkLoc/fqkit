use crate::utils::*;
use bio::io::fastq;
use anyhow::{Result, Ok};
use crossbeam::channel::unbounded;
use log::*;
use std::time::Instant;


pub fn filter_fastq(
    read1: &str,
    read2: &str,
    nbase: usize,
    length: usize,
    complexity: u32,
    average_qual: u8,
    phred: u8,
    chunks: usize,
    ncpu: usize,
    failed: &str,
    out1: &str,
    out2: &str,
    compression_level: u32,
    quiet: bool,
) -> Result<()> {
    let start = Instant::now();
    if !quiet {
        info!("read forward reads from file: {}", read1);
        info!("read reverse reads from file: {}", read2);
        if ![33u8, 64u8].contains(&phred) {
            error!("invalid phred value");
            std::process::exit(1);
        }
        if ncpu <= 1 {
            info!("thread num is: {}", ncpu);
        } else {
            info!("additional thread num is: {}", ncpu);
        }
        info!("output clean read1 file: {}", out1);
        info!("output clean read2 file: {}", out2);
    }

    let fq_reader1 = file_reader(&Some(read1)).map(fastq::Reader::new)?;
    let fq_reader2 = file_reader(&Some(read2)).map(fastq::Reader::new)?;
    let mut out_writer1 = file_writer(&Some(out1), compression_level).map(fastq::Writer::new)?;
    let mut out_writer2 = file_writer(&Some(out2), compression_level).map(fastq::Writer::new)?;
    let mut failed_writer = file_writer(&Some(failed), compression_level).map(fastq::Writer::new)?;
    let complex = complexity as usize;
    let (mut pe_ok, mut pe_fail) = (0usize,0usize);
    let (mut fail_n, mut fail_qual,mut fail_len, mut fail_complex) = (0usize,0usize,0usize,0usize);

    if ncpu <= 1 {
        for (rec1,rec2) in fq_reader1.records().flatten().zip(fq_reader2.records().flatten()) {
            if rec1.seq().iter().filter(|v| v == &&b'N').count() > nbase || rec2.seq().iter().filter(|v| v == &&b'N').count() > nbase {
                failed_writer.write_record(&rec1)?;
                failed_writer.write_record(&rec2)?;
                pe_fail += 1;
                fail_n += 1;
                continue;
            }
            if rec1.seq().len() < length || rec2.seq().len() < length {
                failed_writer.write_record(&rec1)?;
                failed_writer.write_record(&rec2)?;
                pe_fail += 1;
                fail_len += 1;
                continue;
            }
            let complx1 = (rec1.seq().iter().skip(1)
                .zip(rec1.seq().iter())
                .filter(|(q1,q2)| {q1 != q2})
                .count()  as f64 / rec1.seq().len() as f64 * 100.0 ) as usize;
            let complx2 = (rec2.seq().iter().skip(1)
                .zip(rec1.seq().iter())
                .filter(|(q1,q2)| {q1 != q2})
                .count()  as f64 / rec2.seq().len() as f64 * 100.0 ) as usize;
            if complx1 < complex || complx2 < complex {
                failed_writer.write_record(&rec1)?;
                failed_writer.write_record(&rec2)?;
                pe_fail += 1;
                fail_complex += 1;
                continue;
            }
            if phred_mean(rec1.qual(),phred) < average_qual || phred_mean(rec2.qual(), phred) < average_qual {
                failed_writer.write_record(&rec1)?;
                failed_writer.write_record(&rec2)?;
                pe_fail += 1;
                fail_qual += 1;
                continue;
            }
            pe_ok += 1;
            out_writer1.write_record(&rec1)?;
            out_writer2.write_record(&rec2)?;
        }
    } else {
        let mut chunk = chunks;
        if chunk == 0 {
            warn!("pe read conut in chunk can't be: {}, changed to default value.",chunk);
            chunk = 5000;
        }

        let (tx,rx) = unbounded();
        let mut fq_iter1 = fq_reader1.records();
        let mut fq_iter2 = fq_reader2.records();
        loop {
            let pe_vec: Vec<_> = fq_iter1.by_ref().take(chunk).flatten().zip(fq_iter2.by_ref().take(chunk).flatten()).collect();
            if pe_vec.is_empty() {
                break;
            }
            tx.send(pe_vec).unwrap();
        }
        drop(tx);

        crossbeam::scope(|s| {
            let (tx2,rx2) = unbounded();
            let _handles: Vec<_> = (0..ncpu).map(|_| {
                let tx_tmp = tx2.clone();
                let rx_tmp = rx.clone();
                s.spawn(move |_| {
                    for vec_pe in rx_tmp.iter() {
                        let mut passed = vec![];
                        let mut failed = vec![];
                        let (mut fail_n, mut fail_qual,mut fail_len, mut fail_complex) = (0usize,0usize,0usize,0usize);
                        for (rec1,rec2) in vec_pe {
                            if rec1.seq().iter().filter(|v| v == &&b'N').count() > nbase || rec2.seq().iter().filter(|v| v == &&b'N').count() > nbase {
                                failed.push((rec1,rec2));
                                fail_n += 1;
                                continue;
                            }
                            if rec1.seq().len() < length || rec2.seq().len() < length {
                                failed.push((rec1,rec2));
                                fail_len += 1;
                                continue;
                            }
                            let complx1 = (rec1.seq().iter().skip(1)
                                .zip(rec1.seq().iter())
                                .filter(|(q1,q2)| {q1 != q2})
                                .count()  as f64 / rec1.seq().len() as f64 * 100.0 ) as usize;
                            let complx2 = (rec2.seq().iter().skip(1)
                                .zip(rec1.seq().iter())
                                .filter(|(q1,q2)| {q1 != q2})
                                .count()  as f64 / rec2.seq().len() as f64 * 100.0 ) as usize;
                            if complx1 < complex || complx2 < complex {
                                failed.push((rec1,rec2));
                                fail_complex += 1;
                                continue;
                            }
                            if phred_mean(rec1.qual(),phred) < average_qual || phred_mean(rec2.qual(), phred) < average_qual {
                                failed.push((rec1,rec2));
                                fail_qual += 1;
                                continue;
                            }
                            passed.push((rec1,rec2));
                        }
                        tx_tmp.send((passed,failed,fail_n,fail_len,fail_qual,fail_complex)).unwrap();
                    }
                });
            }).collect();
            drop(tx2);

            for (vec_pass,vec_failed,n,len,qual,complex) in rx2.iter() {
                fail_n += n;
                fail_len += len;
                fail_qual += qual;
                fail_complex += complex;

                for (rec1,rec2) in vec_pass {
                    pe_ok += 1;
                    out_writer1.write_record(&rec1).unwrap();
                    out_writer2.write_record(&rec2).unwrap();
                }
                for (rec1,rec2) in vec_failed.iter() {
                    pe_fail += 1;
                    failed_writer.write_record(&rec1).unwrap();
                    failed_writer.write_record(&rec2).unwrap();
                }
            }
        }).unwrap();
    }
    out_writer1.flush()?;
    out_writer2.flush()?;
    failed_writer.flush()?;
    
    if !quiet {
        info!("total clean pe reads number (r1+r2): {}", pe_ok*2);
        info!("total failed pe reads number (r1+r2): {}", pe_fail*2);
        info!("read1 or read2 failed due to low complexity (r1+r2): {}",fail_complex*2);
        info!("read1 or read2 failed due to low quality (r1+r2): {}",fail_n*2);
        info!("read1 or read2 failed due to too many N (r1+r2): {}",fail_n*2);
        info!("read1 or read2 failed due to too short (r1+r2): {}",fail_len*2);
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}


fn phred_mean(
    qual: &[u8],
    phred: u8,
) -> u8 {
    let ave_error = qual
        .iter()
        .map(|x| { 10.0f64.powf((x - phred) as f64 / -10.0) })
        .sum::<f64>() / qual.len() as f64;
    
    (-10.0f64 * ave_error.log10()).round() as u8
}