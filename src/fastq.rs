use bio::io::fastq;
use colored::*;
use rand::{prelude::*, Rng};
use rand_pcg::Pcg64;
use std::{collections::HashMap, io::Result};
use std::io::BufRead;

use crate::utils::*;

// reduce much memory but cost more time
pub fn select_fastq(file: &Option<&str>, n: usize, seed: u64, out: &Option<&str>) -> Result<()> {
    let mut rng = Pcg64::seed_from_u64(seed);
    let mut get: Vec<usize> = Vec::with_capacity(n);

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    for (order, _) in fq_reader.records().flatten().enumerate() {
        if order < n {
            get.push(order);
        } else {
            let ret = rng.gen_range(0..=order);
            if ret < n {
                get[ret] = order;
            }
        }
    }

    let fo = file_writer(out)?;
    let mut w = fastq::Writer::new(fo);

    let fq_reader2 = fastq::Reader::new(file_reader(file)?);
    for (order, rec) in fq_reader2.records().flatten().enumerate() {
        if get.contains(&order) {
            w.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
        }
    }
    Ok(())
}

// fast mode but cost more memory
pub fn select_fastq2(file: &Option<&str>, n: usize, seed: u64, out: &Option<&str>) -> Result<()> {
    let mut rng = Pcg64::seed_from_u64(seed);
    let mut get: Vec<fastq::Record> = Vec::with_capacity(n);

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    for (order, rec) in fq_reader.records().flatten().enumerate() {
        if order < n {
            get.push(rec);
        } else {
            let ret = rng.gen_range(0..=order);
            if ret < n {
                get[ret] = rec;
            }
        }
    }

    let fo = file_writer(out)?;
    let mut w = fastq::Writer::new(fo);
    for rec in get {
        w.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
    }
    Ok(())
}

pub fn fq2fa(file: &Option<&str>, out: &Option<&str>) -> Result<()> {
    let fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut fo = file_writer(out)?;

    for rec in fq_reader.records().flatten() {
        let pre = rec.id();
        let seq = std::str::from_utf8(rec.seq()).expect("Invalid UTF-8 sequence");
        let fa = match rec.desc() {
            Some(desc) => {
                format!(">{} {}\n{}\n", pre, desc, seq)
            }
            None => {
                format!(">{}\n{}\n", pre, seq)
            }
        };
        write!(&mut fo, "{}", fa)?;
    }
    Ok(())
}

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

pub fn stat_fq(inp: &Option<&str>, pre_sum: &str, pre_cyc: &Option<&str>, phred: u8) -> Result<()> {
    if ![33u8, 64u8].contains(&phred) {
        eprintln!("{}", "[error]: invalid phred value".red());
        std::process::exit(1);
    }
    let fq = fastq::Reader::new(file_reader(inp)?);
    let mut fo = file_writer(&Some(pre_sum))?;
    let mut fc = file_writer(pre_cyc)?;

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
                    eprintln!("{}", "[error]: invalid base".red());
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
                        eprintln!("{}", "[error]: invalid base".red());
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
    writeln!(&mut fc, "{}", header.join("\t"))?;

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
        writeln!(&mut fc, "{}", out.join("\t"))?;
    }

    Ok(())
}

pub fn remove_read(
    file: &Option<&str>,
    out: &Option<&str>,
    name: &str,
) -> Result<()> {
    let mut ids = vec![];
    let mut cot = 0usize;
    let list = file_reader(&Some(name))?;
    for i in list.lines().flatten(){
        ids.push(i);
        cot += 1;
    }
    if cot == 0 {
        eprintln!("{}", "[error]: read name list is empty.".red());
        std::process::exit(1);
    }

    let fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut fq_writer = fastq::Writer::new(file_writer(out)?);
    for rec in fq_reader.records().flatten() {
        if !ids.contains(&rec.id().to_string()) {
            fq_writer.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;    
        }    
    }
    Ok(())
}

pub fn split_interleaved(
    file: &Option<&str>,
    out_dir: &str,
    out_pre: &str,
) -> Result<()> {
    let pre1 = format!("{}/{}_r1.fq.gz", out_dir, out_pre);
    let pre2 = format!("{}/{}_r2.fq.gz", out_dir, out_pre);
    let mut fh1 = fastq::Writer::new(file_writer_append(&pre1)?);
    let mut fh2 = fastq::Writer::new(file_writer_append(&pre2)?);
    
    let mut n = 0;
    let fq_reader = fastq::Reader::new(file_reader(file)?);
    for rec in fq_reader.records().flatten() {
        n += 1;
        if n == 1 {
            fh1.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
        } else {
            n = 0;
            fh2.write(rec.id(), rec.desc(), rec.seq(), rec.qual())?;
        }
    }
    Ok(())
}

