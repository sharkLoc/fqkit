use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use nthash::nthash;
use std::collections::HashMap;

pub fn kmer_count(
    input: Option<&String>,
    kmer_len: usize,
    header: bool,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let reader = file_reader(input).map(fastq::Reader::new)?;

    let mut writer = file_writer(output, compression_level, stdout_type)?;
    let mut kmers = HashMap::new();

    for rec in reader.records().map_while(Result::ok) {
        let (mut sidx, mut eidx) = (0, kmer_len);
        let khash = nthash(rec.seq(), kmer_len);
        let len = rec.seq().len();

        while eidx <= len {
            let kseq = &rec.seq()[sidx..eidx];
            let khash_this = nthash(kseq, kmer_len)[0];
            if khash.contains(&khash_this) {
                *kmers.entry(kseq.to_owned()).or_insert(0_u64) += 1;
            }
            sidx += 1;
            eidx += 1;
        }
    }

    if header {
        writer.write_all("kmer\tcount\n".as_bytes())?;
    }
    for (k, v) in kmers {
        writer.write_all(k.as_slice())?;
        writer.write_all(format!("\t{}\n", v).as_bytes())?;
    }
    writer.flush()?;
    Ok(())
}
