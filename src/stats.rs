use bio::io::fastq;
use std::{collections::HashMap, io::Result};
use log::*;
use crate::utils::*;
use std::time::Instant;

struct Nt<'a> {
    seq: &'a [u8],
    qual: &'a [u8],
}
type Cnt = (usize, usize, usize, usize, usize);

impl<'a> Nt<'a> {
    #[inline]
    fn new() -> Self {
        Nt {
            seq: &[],
            qual: &[],
        }
    }
    #[inline]
    fn count_nt(&self) -> Cnt {
        let mut atgcn: Cnt = (0, 0, 0, 0, 0);
        for i in self.seq.iter() {
            if *i == b'A' {
                atgcn.0 += 1;
            }
            if *i == b'T' {
                atgcn.1 += 1;
            }
            if *i == b'G' {
                atgcn.2 += 1;
            }
            if *i == b'C' {
                atgcn.3 += 1;
            }
            if *i == b'N' {
                atgcn.4 += 1;
            }
        }
        atgcn
    }
    #[inline]
    fn q_value(&self, phred: u8) -> (usize, usize) {
        let mut qva = (0usize, 0usize);
        for i in self.qual.iter() {
            if *i - phred >= 20 {
                qva.0 += 1;
            }
            if *i - phred >= 30 {
                qva.1 += 1;
            }
        }
        qva
    }
}

#[allow(non_camel_case_types)]
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
    num_q20: usize,
    rate_q20: f64,
    num_q30: usize,
    rate_q30: f64,
    rate_gc: f64,
    ave_len: f64,
    max_len: usize,
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
            num_q20: 0,
            rate_q20: 0.0,
            num_q30: 0,
            rate_q30: 0.0,
            rate_gc: 0.0,
            ave_len: 0.0,
            max_len: 0,
        }
    }
    fn calc(&mut self) {
        self.ave_len = self.num_base as f64 / self.num_read as f64;
        self.rate_a = self.num_a as f64 / self.num_base as f64;
        self.rate_t = self.num_t as f64 / self.num_base as f64;
        self.rate_g = self.num_g as f64 / self.num_base as f64;
        self.rate_c = self.num_c as f64 / self.num_base as f64;
        self.rate_n = self.num_n as f64 / self.num_base as f64;
        self.rate_gc = (self.num_g + self.num_c) as f64 / self.num_base as f64;
        self.rate_q20 = self.num_q20 as f64 / self.num_base as f64;
        self.rate_q30 = self.num_q30 as f64 / self.num_base as f64;
    }
}

pub fn stat_fq(
    inp: &Option<&str>, 
    pre_sum: &str, 
    pre_cyc: &Option<&str>, 
    phred: u8,
    quiet: bool,
) -> Result<()> {
    if ![33u8, 64u8].contains(&phred) {
        error!("invalid phred value");
        std::process::exit(1);
    }
    let start = Instant::now();
    if !quiet {
        if let Some(inp) = inp {
            info!("reading from file: {}", inp);
        } else {
            info!("reading from stdin");
        }
    }

    let fq = fastq::Reader::new(file_reader(inp)?);
    let mut fo = file_writer(&Some(pre_sum))?;
    let mut fc = file_writer(pre_cyc)?;
    if !quiet {
        info!("summary result write to file: {}",pre_sum);
        if pre_cyc.is_some() {
            info!("cycle result write to file: {}",pre_cyc.unwrap());
        }
    }

    let mut cnt = info::new();
    let mut hs_pos: HashMap<usize, HashMap<i32, usize>> = HashMap::new();
    let mut max_qva = 0;
    for rec in fq.records().flatten() {
        let mut info = Nt::new();
        info.qual = rec.qual();
        info.seq = rec.seq();
        // each cycle info
        for (pos, (sf, sq)) in rec.seq().iter().zip(rec.qual().iter()).enumerate() {
            let pos1 = pos + 1;
            let q = (sq - phred) as i32;
            if q > max_qva {
                max_qva = q;
            }
            let mut hs_ct: HashMap<i32, usize> = HashMap::new();
            match *sf {
                b'A' => {
                    *hs_ct.entry(-5).or_insert(0) += 1;
                }
                b'T' => {
                    *hs_ct.entry(-4).or_insert(0) += 1;
                }
                b'G' => {
                    *hs_ct.entry(-3).or_insert(0) += 1;
                }
                b'C' => {
                    *hs_ct.entry(-2).or_insert(0) += 1;
                }
                b'N' => {
                    *hs_ct.entry(-1).or_insert(0) += 1;
                }
                _ => {
                    error!("invalid base");
                    panic!();
                }
            }
            *hs_ct.entry(q).or_insert(0) += 1;

            if hs_pos.contains_key(&pos1) {
                let exists = hs_pos.get_mut(&pos1).unwrap();
                match *sf {
                    b'A' => {
                        *exists.entry(-5).or_insert(0) += 1;
                    }
                    b'T' => {
                        *exists.entry(-4).or_insert(0) += 1;
                    }
                    b'G' => {
                        *exists.entry(-3).or_insert(0) += 1;
                    }
                    b'C' => {
                        *exists.entry(-2).or_insert(0) += 1;
                    }
                    b'N' => {
                        *exists.entry(-1).or_insert(0) += 1;
                    }
                    _ => {
                        error!("invalid base");
                        panic!();
                    }
                }
                *exists.entry(q).or_insert(0) += 1;
            } else {
                hs_pos.entry(pos1).or_insert(hs_ct);
            }
        }

        let nt: Cnt = info.count_nt();
        let qv = info.q_value(phred);

        cnt.num_read += 1;
        let rlen = rec.seq().len();
        cnt.num_base += rlen;

        if rlen > cnt.max_len {
            cnt.max_len = rlen;
        }

        cnt.num_a += nt.0;
        cnt.num_t += nt.1;
        cnt.num_g += nt.2;
        cnt.num_c += nt.3;
        cnt.num_n += nt.4;
        cnt.num_q20 += qv.0;
        cnt.num_q30 += qv.1;
    }
    cnt.calc();

    // output summary result
    writeln!(&mut fo, "read average length:\t{}", cnt.ave_len)?;
    writeln!(&mut fo, "read max length:\t{}", cnt.max_len)?;
    writeln!(&mut fo, "total gc content(%):\t{:.2}", cnt.rate_gc * 100.0)?;
    writeln!(&mut fo, "total read count:\t{}", cnt.num_read)?;
    writeln!(&mut fo, "total base count:\t{}\n", cnt.num_base)?;
    writeln!(
        &mut fo,
        "base A count:\t{}\t({:.2}%)",
        cnt.num_a,
        cnt.rate_a * 100.0
    )?;
    writeln!(
        &mut fo,
        "base T count:\t{}\t({:.2}%)",
        cnt.num_t,
        cnt.rate_t * 100.0
    )?;
    writeln!(
        &mut fo,
        "base G count:\t{}\t({:.2}%)",
        cnt.num_g,
        cnt.rate_g * 100.0
    )?;
    writeln!(
        &mut fo,
        "base C count:\t{}\t({:.2}%)",
        cnt.num_c,
        cnt.rate_c * 100.0
    )?;
    writeln!(
        &mut fo,
        "base N count:\t{}\t({:.2}%)\n",
        cnt.num_n,
        cnt.rate_n * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 20 or higher (Q20+) (%)\t{}\t({:.2}%)",
        cnt.num_q20,
        cnt.rate_q20 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 30 or higher (Q30+) (%)\t{}\t({:.2}%)",
        cnt.num_q30,
        cnt.rate_q30 * 100.0
    )?;
    fo.flush()?;

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
    //writeln!(&mut fc, "{}", header.join("\t"))?;
    fc.write(header.join("\t").as_bytes())?;
    fc.write(b"\n")?;

    for i in 0..cnt.max_len {
        let idx = i + 1;
        let this = hs_pos.get(&idx).unwrap();
        let mut out = Vec::new();
        out.push(format!("cyc{}", idx));
        for x in -5..=max_qva {
            let n = this.get(&x).unwrap_or(&0);
            if x < 0 {
                let rate = *n as f64
                    / (this.get(&-5).unwrap_or(&0)
                        + this.get(&-4).unwrap_or(&0)
                        + this.get(&-3).unwrap_or(&0)
                        + this.get(&-2).unwrap_or(&0)
                        + this.get(&-1).unwrap_or(&0)) as f64
                    * 100.0;
                out.push(format!("{}:({:.2}%)", n, rate));
            } else {
                out.push(format!("{}", n));
            }
        }
        //writeln!(&mut fc, "{}", out.join("\t"))?;
        fc.write(out.join("\t").as_bytes())?;
        fc.write(b"\n")?;
    }
    fc.flush()?;

    if !quiet {
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}
