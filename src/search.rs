use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use regex::Regex;
use crossbeam::channel::unbounded;
use crate::utils::*;
use std::time::Instant;


pub fn search_fq(
    fq: &Option<&str>,
    pat: &str,
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
        info!("additional thread num is: {}", ncpu);
        if let Some(out) = out {
            info!("reads write to file: {}", out);
        } else {
            info!("reads write to stdout");
        }
    }

    let mut num = 0usize;
    let mut fo = file_writer(out)
        .map(fastq::Writer::new)?;
    let fq_reader = file_reader(fq)
            .map(fastq::Reader::new).unwrap();
    
    if ncpu == 1 || ncpu == 0 {
        let re = Regex::new(pat).unwrap();
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
        for rec in fq_reader.records().flatten() {
            tx.send(rec).unwrap();
        }
        drop(tx);
        crossbeam::scope(|s| {
            /*s.spawn( |_| {
                let fq_reader = file_reader(fq).map(fastq::Reader::new).unwrap();
                for rec in fq_reader.records().flatten() {
                    tx.send(rec).unwrap();
                }
                drop(tx);
            });*/
            let ( tx2, rx2) = unbounded(); //unbounded::<bio::io::fastq::Record>();
            let _handles: Vec<_> = (0..ncpu).map(|_| {
                let tx_tmp = tx2.clone();
                let rx_tmp = rx.clone();
                let re = Regex::new(pat).unwrap();
                s.spawn(move  |_| {    
                    for rec in rx_tmp.iter() {
                        let fq_str = std::str::from_utf8(rec.seq()).unwrap();
                        if  re.is_match(fq_str) {
                            tx_tmp.send(rec).unwrap();
                        }
                    }
                });  
            }).collect();
            drop(tx2);

            for rec in rx2.iter() {
                num += 1;
                fo.write_record(&rec).unwrap();
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