use crate::utils::*;
use anyhow::{Result,Ok};
use bio::io::{fastq, fastq::Record};
use bio::alignment::{pairwise::*, AlignmentOperation};
use crossbeam::channel::unbounded;
use log::*;
use std::time::Instant;


pub fn join_fastq(
    fq1: &str,
    fq2: &str,
    min_overlap: usize,
    max_mismatch: usize,
    chunk: usize,
    ncpu: usize,
    merge_fq: &str,
    out_r1: &str,
    out_r2: &str,
    compression_level: u32,
    quiet: bool,
) -> Result<()> {
    if !quiet {
        info!("read forward reads from file: {}", fq1);
        info!("read reverse reads from file: {}", fq2);
        if ncpu <= 1 {
            info!("thread num is: {}", ncpu);
        } else {
            info!("additional thread num is: {}", ncpu);
        }
        info!("outout unmerged read1 file: {}", out_r1);
        info!("outout unmerged read2 file: {}", out_r2);
        info!("output merged pe reads file: {}",merge_fq);
    }
    let start = Instant::now();

    let fq_reader1 = file_reader(&Some(fq1)).map(fastq::Reader::new)?;
    let fq_reader2 = file_reader(&Some(fq2)).map(fastq::Reader::new)?;
    let mut out_writer1 = file_writer(&Some(out_r1), compression_level).map(fastq::Writer::new)?;
    let mut out_writer2 = file_writer(&Some(out_r2), compression_level).map(fastq::Writer::new)?;
    let mut merge_writer = file_writer(&Some(merge_fq), compression_level).map(fastq::Writer::new)?;
    
    let (mut total_count, mut merge_count, mut min_overlap_fail, mut max_mismatch_fail ) = (0usize, 0usize, 0usize, 0usize);
    let scoring = Scoring {
        gap_open: -5,
        gap_extend: -1,
        match_fn: |a: u8, b: u8| if a == b { 1i32 } else { -1i32 },
        match_scores: Some((1, -1)),
        xclip_prefix: 0,
        xclip_suffix: MIN_SCORE,
        yclip_prefix: 0,
        yclip_suffix: 0,
    };

    if ncpu <= 1 {
        for (rec1,rec2) in fq_reader1.records().flatten().zip(fq_reader2.records().flatten()) {
            total_count += 1;
            let mut aligner = Aligner::with_capacity_and_scoring(rec1.seq().len(), rec2.seq().len(), scoring);
            let rc2_vec = rec2.seq().iter().rev().map(|v| match v {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                _ => b'N',
            }).collect::<Vec<u8>>();
            let rc2 = rc2_vec.as_slice();
            let rc2_qual_vec = rec2.qual().iter().rev().copied().collect::<Vec<u8>>();
            let rc2_qual = rc2_qual_vec.as_slice();
            let alignment = aligner.custom(rec1.seq(), rc2);
            //assert_eq!(rec1.seq().len()-alignment.xstart,alignment.x_aln_len());
            
            let y_clip = alignment.yend - alignment.y_aln_len();
            if min_overlap > alignment.y_aln_len() + y_clip {
                min_overlap_fail += 1;
                if !quiet {
                    trace!("pe reads overlap length {} < min_overlap length: {}",alignment.y_aln_len() + y_clip ,min_overlap)
                }
                out_writer1.write_record(&rec1)?;
                out_writer2.write_record(&rec2)?;
                continue;
            }
            if alignment.xstart > 0  && (alignment.ystart + alignment.y_aln_len() < alignment.ylen) {
                let r1_prefix = &rec1.seq()[0..alignment.xstart];
                let q1_prefix = &rec1.qual()[0..alignment.xstart];
                let r2_suffix = &rc2[alignment.yend..alignment.ylen];
                let q2_suffix = &rc2_qual[alignment.yend..alignment.ylen];
    
                // mismatch in overlap region
                let mut subst = y_clip;
                for each in alignment.operations.iter() {
                    match each {
                        AlignmentOperation::Subst => subst += 1,
                        _ => continue
                    }
                }
                if subst > max_mismatch {
                    max_mismatch_fail += 1;
                    if !quiet { 
                        trace!("mismatch count in overlap region: {} < allowed max mismatch number: {}",subst, max_mismatch); 
                    }
                    out_writer1.write_record(&rec1)?;
                    out_writer2.write_record(&rec2)?;
                    continue;
                }
                // select better quality base for subst in overlap
                merge_count += 1;
                let overlap_r1 = &rec1.seq()[alignment.xstart-y_clip..alignment.xend]; 
                let overlap_q1 = &rec1.qual()[alignment.xstart-y_clip..alignment.xend];
                let overlap_r2 = &rc2[0..alignment.yend];
                let overlap_q2 = &rc2_qual[0..alignment.yend];
                let mut overlap_new_seq = vec![];
                let mut overlap_new_qual = vec![];
                for i in 0..alignment.yend {
                    if overlap_q1[i] >= overlap_q2[i] {
                        overlap_new_seq.push(overlap_r1[i]);
                        overlap_new_qual.push(overlap_q1[i]);
                    } else {
                        overlap_new_qual.push(overlap_q2[i]);
                        overlap_new_seq.push(overlap_r2[i]);
                    }
                }
                let full_seq = [r1_prefix,overlap_new_seq.as_slice(),r2_suffix].concat();
                let full_qual = [q1_prefix,overlap_new_qual.as_slice(),q2_suffix].concat();
                let merged_read = Record::with_attrs(rec1.id(), None, full_seq.as_slice(), full_qual.as_slice());
                merge_writer.write_record(&merged_read)?;
            } 
        }
    } else {
        let mut chunk = chunk;
        if chunk == 0 {
            warn!("pe read conut in chunk can't be: {}, changed to default value.",chunk);
            chunk = 5000;
        }
        let (tx,rx) = unbounded();
        let mut fq_iter1 = fq_reader1.records();
        let mut fq_iter2 = fq_reader2.records();
        loop {
            let pe_vec: Vec<_> = fq_iter1.by_ref().flatten().zip(fq_iter2.by_ref().flatten()).take(chunk).collect(); 
            //fq_iter1.take(chunk).flatten().zip(fq_iter2.take(chunk).flatten()).collect();
            if pe_vec.is_empty() {
                break;
            }
            tx.send(pe_vec).unwrap();
        }
        drop(tx);
        
        crossbeam::scope(|s| {
            let (tx2, rx2) = unbounded();
            let _handles: Vec<_> = (0..ncpu).map(|_| {
                let rx_tmp  = rx.clone();
                let tx_tmp = tx2.clone();
                s.spawn(move |_| {
                    for pe_chunk in rx_tmp {
                        let mut merge_ok = vec![];
                        let mut merge_failed = vec![];
                        let (mut total_count, mut merge_count, mut min_overlap_fail, mut max_mismatch_fail ) = (0usize, 0usize, 0usize, 0usize);
                        for (rec1,rec2) in pe_chunk {
                            total_count += 1;
                            let rc2_vec = rec2.seq().iter().rev().map(|v| match v {
                                b'A' => b'T',
                                b'T' => b'A',
                                b'G' => b'C',
                                b'C' => b'G',
                                _ => b'N',
                            }).collect::<Vec<u8>>();
                            let rc2 = rc2_vec.as_slice();
                            let rc2_qual_vec = rec2.qual().iter().rev().copied().collect::<Vec<u8>>();
                            let rc2_qual = rc2_qual_vec.as_slice();
                            let mut aligner = Aligner::with_capacity_and_scoring(rec1.seq().len(), rec2.seq().len(), scoring);
                            let alignment = aligner.custom(rec1.seq(), rc2);

                            let y_clip = alignment.yend - alignment.y_aln_len();
                            if min_overlap > alignment.y_aln_len() + y_clip {
                                min_overlap_fail += 1;
                                if !quiet {
                                    trace!("pe reads overlap length {} < min_overlap length: {}",alignment.y_aln_len() + y_clip ,min_overlap)
                                }
                                merge_failed.push((rec1,rec2));
                                continue;
                            }
                            if alignment.xstart > 0  && (alignment.ystart + alignment.y_aln_len() < alignment.ylen) {
                                let r1_prefix = &rec1.seq()[0..alignment.xstart];
                                let q1_prefix = &rec1.qual()[0..alignment.xstart];
                                let r2_suffix = &rc2[alignment.yend..alignment.ylen];
                                let q2_suffix = &rc2_qual[alignment.yend..alignment.ylen];
                    
                                // mismatch in overlap region
                                let mut subst = y_clip;
                                for each in alignment.operations.iter() {
                                    match each {
                                        AlignmentOperation::Subst => subst += 1,
                                        _ => continue
                                    }
                                }
                                if subst > max_mismatch {
                                    max_mismatch_fail += 1;
                                    if !quiet { 
                                        trace!("mismatch count in overlap region: {} < allowed max mismatch number: {}",subst, max_mismatch); 
                                    }
                                    merge_failed.push((rec1,rec2));
                                    continue;
                                }
                                // select better quality base for subst in overlap
                                merge_count += 1;
                                let overlap_r1 = &rec1.seq()[alignment.xstart-y_clip..alignment.xend]; 
                                let overlap_q1 = &rec1.qual()[alignment.xstart-y_clip..alignment.xend];
                                let overlap_r2 = &rc2[0..alignment.yend];
                                let overlap_q2 = &rc2_qual[0..alignment.yend];
                                let mut overlap_new_seq = vec![];
                                let mut overlap_new_qual = vec![];
                                for i in 0..alignment.yend {
                                    if overlap_q1[i] >= overlap_q2[i] {
                                        overlap_new_seq.push(overlap_r1[i]);
                                        overlap_new_qual.push(overlap_q1[i]);
                                    } else {
                                        overlap_new_qual.push(overlap_q2[i]);
                                        overlap_new_seq.push(overlap_r2[i]);
                                    }
                                }
                                let full_seq = [r1_prefix,overlap_new_seq.as_slice(),r2_suffix].concat();
                                let full_qual = [q1_prefix,overlap_new_qual.as_slice(),q2_suffix].concat();
                                let merged_read = Record::with_attrs(rec1.id(), None, full_seq.as_slice(), full_qual.as_slice());
                                merge_ok.push(merged_read);
                                //merge_writer.write_record(&merged_read)?;
                            } 
                        }
                        tx_tmp.send((merge_ok,merge_failed,total_count,merge_count, min_overlap_fail, max_mismatch_fail)).unwrap();
                    }    
                });   
            }).collect();
            drop(tx2);

            for (merge_ok,merge_failed,total_count_v,merge_count_v, min_overlap_fail_v, max_mismatch_fail_v) in rx2 {
                for merged in merge_ok {
                    merge_writer.write_record(&merged).unwrap();
                }
                for (rec1,rec2) in merge_failed {
                    out_writer1.write_record(&rec1).unwrap();
                    out_writer2.write_record(&rec2).unwrap();
                }
                merge_count += merge_count_v;
                total_count += total_count_v;
                min_overlap_fail += min_overlap_fail_v;
                max_mismatch_fail += max_mismatch_fail_v;
            }
        }).unwrap();
    }

    merge_writer.flush()?;
    out_writer1.flush()?;
    out_writer2.flush()?;

    if !quiet {
        info!("total pe reads: {}", total_count);
        info!("total merged pe reads: {}",merge_count);
        info!("max_mismatch_failed pe reads: {}",max_mismatch_fail);
        info!("min_overlap_failed pe reads: {}",min_overlap_fail);
        //info!("mean merged read length: {}"); todo
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}