use std::{
    fs::{File, OpenOptions},
    io::{self, prelude::*, BufRead, BufReader, BufWriter, Result, Write},
};

use flate2::{read, write, Compression};

const GZ_MAGIC: [u8; 3] = [0x1f, 0x8b, 0x08];

pub fn is_gzipped(file_name: &str) -> Result<bool> {
    let mut buffer: [u8; 3] = [0; 3];
    let mut fp = File::open(file_name)?;
    let _x = fp.read(&mut buffer)?;
    Ok(buffer == GZ_MAGIC || file_name.ends_with(".gz"))
}

pub fn file_reader(file_in: &Option<&str>) -> Result<Box<dyn BufRead>> {
    if let Some(file_name) = file_in {
        let fp = File::open(file_name)?;
        let flag = is_gzipped(file_name)?;

        if flag {
            Ok(Box::new(BufReader::with_capacity(
                1024 * 1024,
                read::MultiGzDecoder::new(fp),
            )))
        } else {
            Ok(Box::new(BufReader::with_capacity(1024 * 1024, fp)))
        }
    } else {
        let fp = BufReader::new(io::stdin());
        Ok(Box::new(fp))
    }
}

pub fn file_writer(file_out: &Option<&str>) -> Result<Box<dyn Write>> {
    if let Some(file_name) = file_out {
        let fp = File::create(file_name)?;
        if file_name.ends_with(".gz") || file_name.ends_with(".gzip") {
            Ok(Box::new(BufWriter::with_capacity(
                1024 * 1024,
                write::GzEncoder::new(fp, Compression::default()),
            )))
        } else {
            Ok(Box::new(BufWriter::with_capacity(1024 * 1024, fp)))
        }
    } else {
        Ok(Box::new(BufWriter::new(io::stdout())))
    }
}

pub fn file_writer_append(file_out: &str) -> Result<Box<dyn Write>> {
    let fp = OpenOptions::new()
        .append(true)
        .create(true)
        .open(file_out)?;
    
    if file_out.ends_with(".gz") || file_out.ends_with(".gzip") {
        Ok(Box::new(BufWriter::with_capacity(
            1024 * 1024,
            write::GzEncoder::new(fp, Compression::default()),
        )))
    } else {
        Ok(Box::new(BufWriter::with_capacity(1024 * 1024, fp)))    
    }
}
