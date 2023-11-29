use bio::io::fastq;
use std::{collections::HashMap, io::Result, vec};
use log::*;
use crate::utils::*;
use std::time::Instant;


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
        self.num_base = self.num_a + self.num_t + self.num_g+ self.num_c + self.num_n;
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
        info!("summary result write to file: {}",pre_sum);
        if let Some(file) = pre_cyc {
            info!("cycle result write to file: {}", file);
        } else {
            info!("cycle result write to stdout");
        }
    }

    let fq = fastq::Reader::new(file_reader(inp)?);
    let mut fo = file_writer(&Some(pre_sum))?;
    let mut fc = file_writer(pre_cyc)?;

    let mut stat = info::new();
    let mut max_qva = 0;
    let mut min_len: Option<usize>= None;
    let mut each: HashMap<usize, Vec<usize>> = HashMap::new();
    for rec in fq.records().flatten() {

        let this_q = *rec.qual().iter().max().unwrap() - phred;
        if this_q > max_qva { max_qva = this_q; }
        let len = rec.seq().len();
        if len > stat.max_len { stat.max_len = len; }
        
        match min_len {
            Some(v) => { if v>= len { min_len = Some(len) } }
            None => min_len = Some(len) 
        }
     
        let cap = this_q as usize +1 + 5; 
        for (pos, (sf, sq)) in rec.seq().iter().zip(rec.qual().iter()).enumerate() {
            let idx = (sq - phred) as usize;
            if idx >= 20 {
                stat.num_q20 +=1;
                if idx >= 30 { stat.num_q30 +=1; }
            }

            if each.contains_key(&pos) {
                each.get_mut(&pos).unwrap()[idx+5] +=1;
                if sf == &b'A' { each.get_mut(&pos).unwrap()[0] += 1; }
                if sf == &b'T' { each.get_mut(&pos).unwrap()[1] += 1; }
                if sf == &b'G' { each.get_mut(&pos).unwrap()[2] += 1; }
                if sf == &b'C' { each.get_mut(&pos).unwrap()[3] += 1; }
                if sf == &b'N' { each.get_mut(&pos).unwrap()[4] += 1; }
            } else {
                let mut v_tmp = vec![0usize; cap];
                v_tmp[idx+5] += 1;
                if sf == &b'A' { v_tmp[0] += 1; }
                if sf == &b'T' { v_tmp[1] += 1; }
                if sf == &b'G' { v_tmp[2] += 1; }
                if sf == &b'C' { v_tmp[3] += 1; }
                if sf == &b'N' { v_tmp[4] += 1; }
                each.insert(pos, v_tmp.clone());
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
    writeln!(&mut fo, "read average length:\t{}", stat.ave_len)?;
    writeln!(&mut fo, "read min length:\t{}", stat.min_len)?;
    writeln!(&mut fo, "read max length:\t{}", stat.max_len)?;
    writeln!(&mut fo, "total gc content(%):\t{:.2}", stat.rate_gc * 100.0)?;
    writeln!(&mut fo, "total read count:\t{}", stat.num_read)?;
    writeln!(&mut fo, "total base count:\t{}\n", stat.num_base)?;
    writeln!(&mut fo, "base A count:\t{}\t({:.2}%)", stat.num_a, stat.rate_a * 100.0)?;
    writeln!(&mut fo, "base T count:\t{}\t({:.2}%)", stat.num_t, stat.rate_t * 100.0)?;
    writeln!(&mut fo, "base G count:\t{}\t({:.2}%)", stat.num_g, stat.rate_g * 100.0)?;
    writeln!(&mut fo, "base C count:\t{}\t({:.2}%)", stat.num_c, stat.rate_c * 100.0)?;
    writeln!(&mut fo, "base N count:\t{}\t({:.2}%)\n", stat.num_n,stat.rate_n * 100.0)?;
    writeln!(&mut fo, "Number of base calls with quality value of 20 or higher (Q20+) (%)\t{}\t({:.2}%)",stat.num_q20, stat.rate_q20 * 100.0 )?;
    writeln!(&mut fo, "Number of base calls with quality value of 30 or higher (Q30+) (%)\t{}\t({:.2}%)",stat.num_q30, stat.rate_q30 * 100.0 )?;
        
    
    // output cycle result
    let mut header = vec!["cycle".to_string(), "A".to_string(), "T".to_string(), "G".to_string(), "C".to_string(), "N".to_string(),];
    for i in 0..=max_qva {
        header.push(format!("{}", i));
    }
    fc.write(header.join("\t").as_bytes())?;
    fc.write(b"\n")?;

    for x in 0..each.keys().len() { // eq sort cycle 
        let data = each.get(&x).unwrap();
        let mut out = Vec::new();
        out.push(format!("cyc{}", x+1));

        let sum_each = data.iter().take(5).sum::<usize>();
        let index = max_qva as usize + 5;
        for i in 0..=index {
            if i < 5 {
                let rate = data[i] as f64 / sum_each as f64 * 100.0;
                out.push(format!("{}:({:.2}%)", data[i], rate));
            } else {
                if i < index {
                    out.push(format!("{}", data[i]));
                } else {
                    out.push("0".to_string());
                }
            }
        }
        fc.write(out.join("\t").as_bytes())?;
        fc.write(b"\n")?;
    }
    fc.flush()?;

    if !quiet {
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}
