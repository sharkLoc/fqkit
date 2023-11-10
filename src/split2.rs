use crate::utils::*;
use bio::io::fastq;
use anyhow::{Error, Ok};
use log::*;
use std::time::Instant;


pub fn split_chunk(
    file: &Option<&str>,
    num: usize,
    gzip: bool,
    out_pre: &str,
    quiet: bool,
) -> Result<(),Error> {
    let start = Instant::now();
    if !quiet {
        if let Some(file) = file {
            info!("reading from file: {}", file);
        } else {
            info!("reading from stdin");
        }
    }

    let (mut flag, mut index) = (0usize, 0usize);
    let out = if gzip {
        format!("{}{}.fastq.gz",out_pre,index)
    } else {
        format!("{}{}.fastq",out_pre,index)
    };

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut fh = vec![fastq::Writer::new(file_writer(&Some(&out))?)];
    
    if !quiet {
        info!("start to write file: {}",out);
    }
    for rec in fq_reader.records().flatten() {
        if flag < num {
            let fhthis = fh.get_mut(index).unwrap();
            fhthis.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            flag += 1;
        } else {
            index += 1;
            let out = if gzip {
                format!("{}{}.fastq.gz",out_pre,index)
            } else {
                format!("{}{}.fastq",out_pre,index)
            };
            fh.push(fastq::Writer::new(file_writer(&Some(&out))?));
            let fhthis = fh.get_mut(index).unwrap();
            
            if !quiet {
                info!("start to write file: {}",out);
            }
            fhthis.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
            flag = 1; // already write one record in this loop, flag add one  
        }
    }

    if !quiet {
        info!("total chunk number is: {}", index + 1);
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}
