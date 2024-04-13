use crate::utils::*;
use anyhow::{Error, Ok};
use bio::io::fastq;
use log::*;
use std::time::Instant;

pub fn fastq2sam(
    r1: &String,
    r2: Option<&String>,
    sm: &str,
    rg: Option<String>,
    lb: Option<String>,
    pl: Option<String>,
    out: Option<&String>,
    compression_level: u32,
) -> Result<(), Error> {
    let start = Instant::now();
    if let Some(r2) = r2 {
        info!("read file1 from: {}", r1);
        info!("read file2 from: {}", r2);
    } else {
        info!("read file from: {}", r1);
    }
    info!("sample name set: {}", sm);

    let mut sam = file_writer(out, compression_level)?;
    let rg = if let Some(x) = rg {
        x
    } else {
        String::from("A")
    };
    sam.write_fmt(format_args!(
        "@HD\tVN:1.6\tSO:queryname\n@RG\tID:{}\tSM:{}",
        rg, sm
    ))?;
    if let Some(lb) = lb {
        sam.write_fmt(format_args!("\tLB:{}", lb))?;
        info!("library name set: {}", lb);
    }
    if let Some(pl) = pl {
        sam.write_fmt(format_args!("\tPL:{}", pl))?;
        info!("platform set: {}", pl);
    }
    sam.write_all("\n".as_bytes())?;

    if let Some(r2) = r2 {
        let fq1 = file_reader(Some(r1)).map(fastq::Reader::new)?;
        let fq2 = file_reader(Some(r2)).map(fastq::Reader::new)?;
        for (rec1, rec2) in fq1.records().flatten().zip(fq2.records().flatten()) {
            sam.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tRG:Z:{}\n",
                rec1.id(),
                77,
                "*",
                "0",
                "0",
                "*",
                "*",
                "0",
                "0",
                std::str::from_utf8(rec1.seq())?,
                std::str::from_utf8(rec1.qual())?,
                rg,
            ))?;
            sam.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tRG:Z:{}\n",
                rec2.id(),
                141,
                "*",
                "0",
                "0",
                "*",
                "*",
                "0",
                "0",
                std::str::from_utf8(rec2.seq())?,
                std::str::from_utf8(rec2.qual())?,
                rg,
            ))?;
        }
    } else {
        let fq = file_reader(Some(r1)).map(fastq::Reader::new)?;
        for rec in fq.records().flatten() {
            sam.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tRG:Z:{}\n",
                rec.id(),
                4,
                "*",
                "0",
                "0",
                "*",
                "*",
                "0",
                "0",
                unsafe { std::str::from_utf8_unchecked(rec.seq()) },
                unsafe { std::str::from_utf8_unchecked(rec.qual()) },
                rg,
            ))?;
        }
    }
    sam.flush()?;

    info!("time elapsed is: {:?}", start.elapsed());
    Ok(())
}
