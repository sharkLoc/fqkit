use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use crossbeam::channel::bounded;
use log::*;
use regex::RegexBuilder;

#[allow(clippy::too_many_arguments)]
pub fn search_fq(
    fq: Option<&String>,
    pat: &str,
    case: bool,
    invert_match: bool,
    chunk: usize,
    out: Option<&String>,
    ncpu: usize,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {

    let fq_reader = file_reader(fq).map(fastq::Reader::new).unwrap();
    if let Some(file) = fq {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    info!("regex pattern is: {}", pat);

    let mut chunk = chunk;
    if chunk == 0 {
        chunk = 5000;
        warn!(
            "read conut in chunk can't be 0, changed to default value: {}",
            chunk
        );
    }
    let mut num = 0usize;
    let mut fo = file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;
    if let Some(out) = out {
        info!("reads write to file: {}", out);
    } else {
        info!("reads write to stdout");
    }

    if ncpu == 1 || ncpu == 0 {
        let re = RegexBuilder::new(pat)
            .case_insensitive(case)
            .build()
            .unwrap();
        for rec in fq_reader.records().map_while(Result::ok) {
            let fq_str = std::str::from_utf8(rec.seq()).unwrap();
            if invert_match {
                if !re.is_match(fq_str) {
                    num += 1;
                    fo.write_record(&rec)?;
                }
            } else if re.is_match(fq_str) {
                num += 1;
                fo.write_record(&rec)?;
            }
        }
        fo.flush()?;
    } else {
        let (tx, rx) = bounded(5000);
        let mut fqiter = fq_reader.records();
        loop {
            let chunks: Vec<_> = fqiter.by_ref().take(chunk).map_while(Result::ok).collect();
            if chunks.is_empty() {
                break;
            }
            tx.send(chunks).unwrap();
        }
        drop(tx);

        crossbeam::scope(|s| {
            let (tx2, rx2) = bounded(5000);
            let _handles: Vec<_> = (0..ncpu)
                .map(|_| {
                    let tx_tmp = tx2.clone();
                    let rx_tmp = rx.clone();
                    let re = RegexBuilder::new(pat)
                        .case_insensitive(case)
                        .build()
                        .unwrap();

                    s.spawn(move |_| {
                        for vrec in rx_tmp {
                            let mut matchs = vec![];
                            let mut count = 0;
                            for rec in vrec {
                                let fq_str = std::str::from_utf8(rec.seq()).unwrap();
                                if invert_match {
                                    if !re.is_match(fq_str) {
                                        matchs.push(rec);
                                        count += 1;
                                    }
                                } else if re.is_match(fq_str) {
                                    matchs.push(rec);
                                    count += 1;
                                }
                            }
                            tx_tmp.send((matchs, count)).unwrap();
                        }
                    });
                })
                .collect();
            drop(tx2);

            for (recs, count) in rx2.iter() {
                num += count;
                for rec in recs {
                    fo.write_record(&rec).unwrap();
                }
                fo.flush().unwrap();
            }
        })
        .unwrap();
    }

    info!("total reads number: {}", num);
    Ok(())
}
