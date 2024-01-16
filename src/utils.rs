use std::{
    fs::{File, OpenOptions},
    io::{self, prelude::*, BufRead, BufReader, BufWriter, Write},
};
use anyhow::{Ok,Error,Result};


const GZ_MAGIC: [u8; 3] = [0x1f, 0x8b, 0x08];
const BZ_MAGIC: [u8; 3] = [0x42, 0x5a, 0x68];
const XZ_MAGIC: [u8; 6] = [0xfd, 0x37, 0x7a, 0x58, 0x5A, 0x00];
//const ZSTD_MAGIC: [u8; 4] = [0x28, 0xb5, 0x2f, 0xfd]; // TODO
const MAGIC_MAX_LEN: usize  = 6;
const BUFF_SIZE: usize = 1024 * 1024;


fn magic_num(file_name: &str) -> Result<[u8; MAGIC_MAX_LEN], Error> {
    let mut buffer: [u8; MAGIC_MAX_LEN] = [0; MAGIC_MAX_LEN];
    let mut fp = File::open(file_name)?;
    let _ = fp.read(&mut buffer)?;
    Ok(buffer)
}

fn is_gzipped(file_name: &str) -> Result<bool> {
    let buffer = magic_num(file_name)?;
    let gz_or_not = buffer[0] == GZ_MAGIC[0] && buffer[1] == GZ_MAGIC[1] && buffer[2] == GZ_MAGIC[2]; 
    Ok(gz_or_not || file_name.ends_with(".gz"))
}

fn is_bzipped(file_name: &str) -> Result<bool> {
    let buffer = magic_num(file_name)?;
    let bz_or_not = buffer[0] == BZ_MAGIC[0] && buffer[1] == BZ_MAGIC[1] && buffer[2] == BZ_MAGIC[2];
    Ok(bz_or_not || file_name.ends_with(".bz2"))
}

fn is_xz(file_name: &str) -> Result<bool> {
    let buffer = magic_num(file_name)?;
    let xz_or_not = 
           buffer[0] == XZ_MAGIC[0] && buffer[1] == XZ_MAGIC[1] 
        && buffer[2] == XZ_MAGIC[2] && buffer[3] == XZ_MAGIC[3] 
        && buffer[4] == XZ_MAGIC[4] && buffer[5] == XZ_MAGIC[5];
    Ok(xz_or_not || file_name.ends_with(".xz"))
}


pub fn file_reader(file_in: &Option<&str>) -> Result<Box<dyn BufRead>> {
    if let Some(file_name) = file_in {
        let fp = File::open(file_name)?;
        let gz_flag = is_gzipped(file_name)?;
        let bz_flag = is_bzipped(file_name)?;
        let zx_flag = is_xz(file_name)?;

        if gz_flag {
            Ok(Box::new(BufReader::with_capacity(
                BUFF_SIZE,
                flate2::read::MultiGzDecoder::new(fp),
            )))
        } else if bz_flag {
            Ok(Box::new(BufReader::with_capacity(
                BUFF_SIZE, 
                bzip2::read::MultiBzDecoder::new(fp),
            )))
        } else if zx_flag {
            Ok(Box::new(BufReader::with_capacity(
                BUFF_SIZE, 
                xz2::read::XzDecoder::new(fp),
            )))
        } else {
            Ok(Box::new(BufReader::with_capacity(
                BUFF_SIZE, 
                fp
            )))
        }
    } else {
        let fp = BufReader::new(io::stdin());
        Ok(Box::new(fp))
    }
}


pub fn file_writer(
    file_out: &Option<&str>,
    compression_level: u32,
) -> Result<Box<dyn Write>> {
    if let Some(file_name) = file_out {
        let fp = File::create(file_name)?;
        if file_name.ends_with(".gz") {
            Ok(Box::new(BufWriter::with_capacity(
                BUFF_SIZE,
                flate2::write::GzEncoder::new(fp, flate2::Compression::new(compression_level)),
            )))
        } else if file_name.ends_with(".bz2") {
            Ok(Box::new(BufWriter::with_capacity(
                BUFF_SIZE, 
                bzip2::write::BzEncoder::new(fp, bzip2::Compression::new(compression_level))
            )))
        } else if file_name.ends_with(".xz") {
            Ok(Box::new(BufWriter::with_capacity(
                BUFF_SIZE, 
                xz2::write::XzEncoder::new(fp, compression_level)
            )))
        } else {
            Ok(Box::new(BufWriter::with_capacity(BUFF_SIZE, fp)))
        }
    } else {
        Ok(Box::new(BufWriter::new(io::stdout())))
    }
}


pub fn file_writer_append(
    file_out: &str,
    compression_level: u32,
) -> Result<Box<dyn Write>> {
    let fp = OpenOptions::new()
        .append(true)
        .create(true)
        .open(file_out)?;
    
    if file_out.ends_with(".gz") {
        Ok(Box::new(BufWriter::with_capacity(
            BUFF_SIZE,
            flate2::write::GzEncoder::new(fp, flate2::Compression::new(compression_level)),
        )))
    } else if file_out.ends_with(".bz2") {
        Ok(Box::new(BufWriter::with_capacity(
            BUFF_SIZE,
            bzip2::write::BzEncoder::new(fp, bzip2::Compression::new(compression_level)),
        )))
    } else if file_out.ends_with("xz") {
        Ok(Box::new(BufWriter::with_capacity(
            BUFF_SIZE, 
            xz2::write::XzEncoder::new(fp, compression_level)
        )))
    } else {
        Ok(Box::new(BufWriter::with_capacity(
            BUFF_SIZE, 
            fp
        )))    
    }
}
