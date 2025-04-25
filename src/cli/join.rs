use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::info;

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => *b, // unknow base like N
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn join_overlap(
    read1: &str,
    read2: &str,
    max_mismatch_rate: f64,
    min_overlap_len: usize,
    overlap_merge: Option<&String>,
    nonoverlap_pe: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let reader1 = fastq::Reader::new(file_reader(Some(read1))?);
    let reader2 = fastq::Reader::new(file_reader(Some(read2))?);

    let mut count_join = 0usize;
    let mut count_total = 0usize;
    let mut count_base_overlap = 0usize;
    let mut count_miss_overlap = 0usize;
    let mut writer_single =
        fastq::Writer::new(file_writer(overlap_merge, compression_level, stdout_type)?);
    let mut nuoverlap_writer =
        fastq::Writer::new(file_writer(nonoverlap_pe, compression_level, stdout_type)?);

    for (rec1, rec2) in reader1
        .records()
        .map_while(Result::ok)
        .zip(reader2.records().map_while(Result::ok))
    {
        count_total += 1;
        let max_overlap_len = rec1.seq().len().min(rec2.seq().len());
        let rec2_seq_rev = reverse_complement(rec2.seq());

        let mut fine_overlap_len = 0;
        if min_overlap_len <= max_overlap_len {
            for overlap_len in (min_overlap_len..=max_overlap_len).rev() {
                assert!(rec1.seq().len() >= overlap_len && rec2_seq_rev.len() >= overlap_len);
                // overlap region in PE reads
                let over1 = &rec1.seq()[rec1.seq().len() - overlap_len..];
                let over2 = &rec2_seq_rev[..overlap_len];

                let mismatch = over1
                    .iter()
                    .zip(over2.iter())
                    .filter(|(x, y)| x != y)
                    .count();
                let mismatch_rate = mismatch as f64 / overlap_len as f64;

                if max_mismatch_rate < mismatch_rate {
                    continue;
                } // mismatch count too much
                if overlap_len < min_overlap_len {
                    continue;
                } // overlap length is too short
                // pe reads overlaped
                fine_overlap_len = overlap_len;
                count_miss_overlap += mismatch;
                break;
            }
        }

        // build longer single read
        if fine_overlap_len > 0 {
            count_join += 1;
            count_base_overlap += fine_overlap_len;
            let mut single_seq = vec![];
            let mut single_qual = vec![];
            single_seq.extend_from_slice(&rec1.seq()[..rec1.seq().len() - fine_overlap_len]);
            single_qual.extend_from_slice(&rec1.qual()[..rec1.qual().len() - fine_overlap_len]);

            let overlap_r1_qual = &rec1.qual()[rec1.qual().len() - fine_overlap_len..];
            let overlap_r1_seq = &rec1.seq()[rec1.seq().len() - fine_overlap_len..];

            let rec2_qual_rev = rec2.qual().iter().copied().rev().collect::<Vec<u8>>();
            let overlap_r2_qual = &rec2_qual_rev[..fine_overlap_len];
            let overlap_r2_seq = &rec2_seq_rev[..fine_overlap_len];

            for i in 0..fine_overlap_len {
                // Prefer R1 if quality is equal
                if overlap_r1_qual[i] >= overlap_r2_qual[i] {
                    single_seq.push(overlap_r1_seq[i]);
                    single_qual.push(overlap_r1_qual[i]);
                } else {
                    single_seq.push(overlap_r2_seq[i]);
                    single_qual.push(overlap_r2_qual[i]);
                };
            }
            single_seq.extend_from_slice(&rec2_seq_rev[fine_overlap_len..]);
            single_qual.extend_from_slice(&rec2_qual_rev[fine_overlap_len..]);

            let read_long =
                fastq::Record::with_attrs(rec1.id(), rec1.desc(), &single_seq, &single_qual);
            writer_single.write_record(&read_long)?;
        } else {
            // no overlap
            if nonoverlap_pe.is_some() {
                nuoverlap_writer.write_record(&rec1)?;
                nuoverlap_writer.write_record(&rec2)?;
            }
        }
    }
    writer_single.flush()?;
    nuoverlap_writer.flush()?;

    let rate = count_join as f64 / count_total as f64;
    info!("total bases in overlap region: {}", count_base_overlap);
    info!("total mismatchs in overlap region: {}", count_miss_overlap);
    info!("total pe reads overlaped and joined number: {}", count_join);
    info!("total pe reads number: {}", count_total);
    info!("pe reads overlap rate: {:.2}%", rate * 100.0);

    Ok(())
}
