use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::{error, info};
use std::{collections::HashMap, vec};

#[allow(non_camel_case_types)]
#[derive(Debug)]
struct info {
    num_a: usize,
    rate_a: f64,
    num_t: usize,
    rate_t: f64,
    num_g: usize,
    rate_g: f64,
    num_c: usize,
    rate_c: f64,
    num_n: usize,
    rate_n: f64,
    num_read: usize,
    num_base: usize,
    num_q5: usize,
    rate_q5: f64,
    num_q10: usize,
    rate_q10: f64,
    num_q15: usize,
    rate_q15: f64,
    num_q20: usize,
    rate_q20: f64,
    num_q30: usize,
    rate_q30: f64,
    rate_gc: f64,
    ave_len: f64,
    max_len: usize,
    min_len: usize,
}

impl info {
    fn new() -> Self {
        info {
            num_a: 0,
            rate_a: 0.0,
            num_t: 0,
            rate_t: 0.0,
            num_g: 0,
            rate_g: 0.0,
            num_c: 0,
            rate_c: 0.0,
            num_n: 0,
            rate_n: 0.0,
            num_read: 0,
            num_base: 0,
            num_q5: 0,
            rate_q5: 0.0,
            num_q10: 0,
            rate_q10: 0.0,
            num_q15: 0,
            rate_q15: 0.0,
            num_q20: 0,
            rate_q20: 0.0,
            num_q30: 0,
            rate_q30: 0.0,
            rate_gc: 0.0,
            ave_len: 0.0,
            max_len: 0,
            min_len: 0,
        }
    }
    fn calc(&mut self) {
        self.num_base = self.num_a + self.num_t + self.num_g + self.num_c + self.num_n;
        self.ave_len = self.num_base as f64 / self.num_read as f64;
        self.rate_a = self.num_a as f64 / self.num_base as f64;
        self.rate_t = self.num_t as f64 / self.num_base as f64;
        self.rate_g = self.num_g as f64 / self.num_base as f64;
        self.rate_c = self.num_c as f64 / self.num_base as f64;
        self.rate_n = self.num_n as f64 / self.num_base as f64;
        self.rate_gc = (self.num_g + self.num_c) as f64 / self.num_base as f64;
        self.rate_q5 = self.num_q5 as f64 / self.num_base as f64;
        self.rate_q10 = self.num_q10 as f64 / self.num_base as f64;
        self.rate_q15 = self.num_q15 as f64 / self.num_base as f64;
        self.rate_q20 = self.num_q20 as f64 / self.num_base as f64;
        self.rate_q30 = self.num_q30 as f64 / self.num_base as f64;
    }
}

pub fn stat_fq(
    inp: Option<&String>,
    pre_sum: &String,
    pre_cyc: Option<&String>,
    phred: u8,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    if ![33u8, 64u8].contains(&phred) {
        error!("{}", FqkitError::InvalidPhredValue);
        std::process::exit(1);
    }

    let fq = fastq::Reader::new(file_reader(inp)?);
    info!("summary result write to file: {}", pre_sum);
    if let Some(file) = pre_cyc {
        info!("cycle result write to file: {}", file);
    } else {
        info!("cycle result write to stdout");
    }

    let mut fo = file_writer(Some(pre_sum), compression_level, stdout_type)?;
    let mut fc = file_writer(pre_cyc, compression_level, stdout_type)?;

    let mut stat = info::new();
    let mut max_qva = 0;
    let mut min_len: Option<usize> = None;
    let mut each: HashMap<usize, Vec<usize>> = HashMap::new();
    for rec in fq.records().map_while(Result::ok) {
        let this_q = *rec.qual().iter().max().unwrap() - phred;
        if this_q > max_qva {
            max_qva = this_q;
        }
        let len = rec.seq().len();
        if len > stat.max_len {
            stat.max_len = len;
        }

        match min_len {
            Some(v) => {
                if v >= len {
                    min_len = Some(len)
                }
            }
            None => min_len = Some(len),
        }

        for (pos, (sf, sq)) in rec.seq().iter().zip(rec.qual().iter()).enumerate() {
            let idx = (sq - phred) as usize;
            if idx >= 5 {
                stat.num_q5 += 1;
            }
            if idx >= 10 {
                stat.num_q10 += 1;
            }
            if idx >= 15 {
                stat.num_q15 += 1;
            }
            if idx >= 20 {
                stat.num_q20 += 1;
                if idx >= 30 {
                    stat.num_q30 += 1;
                }
            }

            if let std::collections::hash_map::Entry::Vacant(e) = each.entry(pos) {
                let cap = this_q as usize + 1 + 5;
                let mut v_tmp = vec![0usize; cap];
                v_tmp[idx + 5] += 1;
                if sf == &b'A' {
                    v_tmp[0] += 1;
                }
                if sf == &b'T' {
                    v_tmp[1] += 1;
                }
                if sf == &b'G' {
                    v_tmp[2] += 1;
                }
                if sf == &b'C' {
                    v_tmp[3] += 1;
                }
                if sf == &b'N' {
                    v_tmp[4] += 1;
                }
                e.insert(v_tmp);
            } else {
                let tf = each.get_mut(&pos).unwrap();
                let gap = (idx + 5 + 1) as i32 - tf.len() as i32;
                if gap > 0 {
                    for _ in 0..gap {
                        tf.insert(tf.len(), 0);
                    }
                }
                tf[idx + 5] += 1;
                if sf == &b'A' {
                    tf[0] += 1;
                }
                if sf == &b'T' {
                    tf[1] += 1;
                }
                if sf == &b'G' {
                    tf[2] += 1;
                }
                if sf == &b'C' {
                    tf[3] += 1;
                }
                if sf == &b'N' {
                    tf[4] += 1;
                }
            }
        }
    }

    for x in 0..each.len() {
        stat.num_a += each.get(&x).unwrap()[0];
        stat.num_t += each.get(&x).unwrap()[1];
        stat.num_g += each.get(&x).unwrap()[2];
        stat.num_c += each.get(&x).unwrap()[3];
        stat.num_n += each.get(&x).unwrap()[4];
    }
    stat.min_len = min_len.unwrap();
    stat.num_read = each.get(&0).unwrap().iter().take(5).sum::<usize>();
    stat.calc();

    // output summary result
    writeln!(&mut fo, "read average length:\t{:.0}", stat.ave_len)?;
    writeln!(&mut fo, "read min length:\t{}", stat.min_len)?;
    writeln!(&mut fo, "read max length:\t{}", stat.max_len)?;
    writeln!(&mut fo, "total gc content(%):\t{:.2}", stat.rate_gc * 100.0)?;
    writeln!(&mut fo, "total read count:\t{}", stat.num_read)?;
    writeln!(&mut fo, "total base count:\t{}\n", stat.num_base)?;
    writeln!(
        &mut fo,
        "base A count:\t{}\t({:.2}%)",
        stat.num_a,
        stat.rate_a * 100.0
    )?;
    writeln!(
        &mut fo,
        "base T count:\t{}\t({:.2}%)",
        stat.num_t,
        stat.rate_t * 100.0
    )?;
    writeln!(
        &mut fo,
        "base G count:\t{}\t({:.2}%)",
        stat.num_g,
        stat.rate_g * 100.0
    )?;
    writeln!(
        &mut fo,
        "base C count:\t{}\t({:.2}%)",
        stat.num_c,
        stat.rate_c * 100.0
    )?;
    writeln!(
        &mut fo,
        "base N count:\t{}\t({:.2}%)\n",
        stat.num_n,
        stat.rate_n * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 5 or higher (Q5+) (%)\t{}\t({:.2}%)",
        stat.num_q5,
        stat.rate_q5 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 10 or higher (Q10+) (%)\t{}\t({:.2}%)",
        stat.num_q10,
        stat.rate_q10 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 15 or higher (Q15+) (%)\t{}\t({:.2}%)",
        stat.num_q15,
        stat.rate_q15 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 20 or higher (Q20+) (%)\t{}\t({:.2}%)",
        stat.num_q20,
        stat.rate_q20 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 30 or higher (Q30+) (%)\t{}\t({:.2}%)",
        stat.num_q30,
        stat.rate_q30 * 100.0
    )?;

    // output cycle result
    let mut header = vec![
        "cycle".to_string(),
        "A".to_string(),
        "T".to_string(),
        "G".to_string(),
        "C".to_string(),
        "N".to_string(),
    ];
    for i in 0..=max_qva {
        header.push(format!("{}", i));
    }
    fc.write_all(header.join("\t").as_bytes())?;
    fc.write_all(b"\n")?;

    for x in 0..each.keys().len() {
        // eq sort cycle
        let data = each.get(&x).unwrap();
        let mut out = Vec::new();
        out.push(format!("cyc{}", x + 1));

        let sum_each = data.iter().take(5).sum::<usize>();
        let index = max_qva as usize + 5 + 1;
        for i in 0..index {
            if i < 5 {
                let rate = data[i] as f64 / sum_each as f64 * 100.0;
                out.push(format!("{}:({:.2}%)", data[i], rate));
            } else {
                let num = data.get(i).unwrap_or(&0); //or(Some(&0)).unwrap();
                out.push(format!("{}", num));
            }
        }
        fc.write_all(out.join("\t").as_bytes())?;
        fc.write_all(b"\n")?;
    }
    fc.flush()?;

    Ok(())
}
