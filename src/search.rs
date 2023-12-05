use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use regex::RegexBuilder;
use crossbeam::channel::unbounded;
use crate::utils::*;
use std::time::Instant;


pub fn search_fq(
    fq: &Option<&str>,
    pat: &str,
    case: bool,
    chunk: usize,
    out: &Option<&str>,
    ncpu: usize,
    quiet: bool,
) -> Result<(), Error> {
    let start = Instant::now();
    if !quiet {
        if let Some(file) = fq {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
        info!("regex pattern is: {}",pat);
        if ncpu == 1 || ncpu == 0 {
            info!("thread num is: {}", ncpu);
        } else {
            info!("additional thread num is: {}", ncpu);
        }
        if let Some(out) = out {
            info!("reads write to file: {}", out);
        } else {
            info!("reads write to stdout");
        }
    }

    let mut chunk = chunk;
    if chunk == 0 {
        chunk = 5000;
        warn!("read conut in chunk can't be 0, changed to default value: {}",chunk);
    }
    let mut num = 0usize;
    let mut fo = file_writer(out)
        .map(fastq::Writer::new)?;
    let fq_reader = file_reader(fq)
            .map(fastq::Reader::new).unwrap();
    
    if ncpu == 1 || ncpu == 0 {
        let re = RegexBuilder::new(pat).case_insensitive(case).build().unwrap();
        for rec in fq_reader.records().flatten() {
            let fq_str = std::str::from_utf8(rec.seq()).unwrap();
            if  re.is_match(fq_str) {
                num += 1;
                fo.write_record(&rec)?;
            }
        }
        fo.flush()?;
    } else {
        let (tx, rx) = unbounded();
        let mut fqiter = fq_reader.records();
        loop {
            let chunks: Vec<_> = fqiter.by_ref().take(chunk).flatten().collect();
            if chunks.is_empty() { break; }
            tx.send(chunks).unwrap();
        }
        drop(tx);
        
        crossbeam::scope(|s| {
            let ( tx2, rx2) = unbounded(); 
            let _handles: Vec<_> = (0..ncpu).map(|_| {
                let tx_tmp = tx2.clone();
                let rx_tmp = rx.clone();
                let re = RegexBuilder::new(pat).case_insensitive(case).build().unwrap();
                s.spawn(move  |_| {    
                    for vrec in rx_tmp {
                        let mut matchs = vec![];
                        let mut count = 0;
                        for rec in vrec {
                            let fq_str = std::str::from_utf8(rec.seq()).unwrap();
                            if  re.is_match(fq_str) {
                                matchs.push(rec);
                                count += 1;
                            }
                        }
                        tx_tmp.send((matchs, count)).unwrap();
                    }
                });  
            }).collect();
            drop(tx2);

            for (recs, count) in rx2.iter() {
                num += count;
                for rec in recs {
                    fo.write_record(&rec).unwrap();
                }
                fo.flush().unwrap();
            }
        }).unwrap();
    }
    
    if !quiet {
        info!("total reads matched number: {}",num);
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}