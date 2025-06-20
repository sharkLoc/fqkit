use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use paraseq::fastq;
use std::collections::HashMap;

pub fn kmer_count(
    input: Option<&String>,
    kmer_len: usize,
    header: bool,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut reader = file_reader(input).map(fastq::Reader::new)?;
    let mut rset = fastq::RecordSet::default();

    let mut writer = file_writer(output, compression_level, stdout_type)?;
    let mut kmers = HashMap::new();

    while rset.fill(&mut reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            let seq = rec.seq();
            if seq.len() < kmer_len {
                continue;
            }
            for kmer in seq.windows(kmer_len) {
                *kmers.entry(kmer.to_vec()).or_insert(0_u64) += 1;
            }
        }
    }

    if header {
        writer.write_all(b"kmer\tcount\n")?;
    }
    for (k, v) in kmers {
        writer.write_all(&k)?;
        writer.write_all(format!("\t{}\n", v).as_bytes())?;
    }
    writer.flush()?;
    Ok(())
}
