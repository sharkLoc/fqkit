use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::info;
use paraseq::fastq;
use rand::prelude::*;
use rand_pcg::Pcg64;
use std::collections::HashMap;

pub fn shuffle_fastq(
    file: Option<&String>,
    seed: u64,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut rng = Pcg64::seed_from_u64(seed);
    let mut fq_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut rset = fastq::RecordSet::default();
    info!("rand seed: {}", seed);

    let mut reads_map = HashMap::new();
    let mut index = 0usize;

    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            reads_map.insert(
                index,
                vec![
                    rec.id().to_owned(),
                    rec.seq().to_owned(),
                    rec.qual().to_owned(),
                ],
            );
            index += 1;
        }
    }

    info!("all records has been readed into memory, start shuffle ...");
    let mut shuffled_indices: Vec<usize> = (0..index).collect();
    shuffled_indices.shuffle(&mut rng);

    info!("shuffle done, start write to output ...");
    let mut writer = file_writer(out, compression_level, stdout_type)?;
    for idx in shuffled_indices {
        if let Some(reads) = reads_map.get(&idx) {
            write_record(&mut writer, reads[0].as_slice(), &reads[1], &reads[2])?;
        }
    }
    writer.flush()?;

    info!("shuffle completed, output written successfully.");
    Ok(())
}
