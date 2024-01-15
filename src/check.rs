use crate::utils::*;
use bio::io::fastq;
use anyhow::Result;
use log::*;
use std::time::Instant;


pub fn check_fastq(
    file: &Option<&str>,
    save: bool,
    out: &Option<&str>,
    compression_level: u32,
) -> Result<()> {

    if let Some(file) = file {
        info!("reading from file: {}", file);
    } else {
        info!("reading from stdin");
    }

    let start = Instant::now();
    let (mut total,mut ok_read, mut err_read) = (0,0,0);
    let fp_reader = file_reader(file).map(fastq::Reader::new)?;
    if save {
        let mut out_writer = file_writer(out, compression_level).map(fastq::Writer::new)?;
        for rec in fp_reader.records().flatten() {
            total +=1;
            match rec.check() {
                Ok(_) => {
                    ok_read += 1;
                    out_writer.write_record(&rec)?;
                }
                Err(s) => {
                    err_read += 1;
                    error!("record order: {}, id: {}, {}", total,rec.id(), s);
                }
            };
        }
        out_writer.flush()?;
    } else {
        for rec in fp_reader.records().flatten() {
            total += 1;
            match rec.check() {
                Ok(_) => {
                    ok_read += 1;
                }
                Err(s) => {
                    err_read += 1;
                    error!("record order: {}, id: {}, {}", total,rec.id(), s);
                }
            };
        }
    }


    info!("total reads num: {}, ok reads number: {}, error reads number: {}", total, ok_read, err_read);
    info!("time elapsed is: {:?}",start.elapsed()); 
    
    Ok(())
}
