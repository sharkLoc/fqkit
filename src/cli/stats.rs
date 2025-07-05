use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{error, info};
use paraseq::{
    fastq,
    fastx::Record,
    parallel::{ParallelProcessor, ParallelReader, ProcessError},
};
use parking_lot::Mutex;
use std::sync::Arc;
use std::{collections::HashMap, vec};

#[allow(non_camel_case_types)]
#[derive(Debug, Clone)]
struct Info {
    num_a: usize,
    total_num_a: Arc<Mutex<usize>>,
    rate_a: f64,
    num_t: usize,
    total_num_t: Arc<Mutex<usize>>,
    rate_t: f64,
    num_g: usize,
    total_num_g: Arc<Mutex<usize>>,
    rate_g: f64,
    num_c: usize,
    total_num_c: Arc<Mutex<usize>>,
    rate_c: f64,
    num_n: usize,
    total_num_n: Arc<Mutex<usize>>,
    rate_n: f64,
    num_read: usize,
    total_num_read: Arc<Mutex<usize>>,
    num_base: usize,
    total_num_base: Arc<Mutex<usize>>,
    num_q5: usize,
    total_num_q5: Arc<Mutex<usize>>,
    rate_q5: f64,
    num_q10: usize,
    total_num_q10: Arc<Mutex<usize>>,
    rate_q10: f64,
    num_q15: usize,
    total_num_q15: Arc<Mutex<usize>>,
    rate_q15: f64,
    num_q20: usize,
    total_num_q20: Arc<Mutex<usize>>,
    rate_q20: f64,
    num_q30: usize,
    total_num_q30: Arc<Mutex<usize>>,
    rate_q30: f64,
    rate_gc: f64,
    ave_len: f64,
    max_len: usize,
    total_max_len: Arc<Mutex<usize>>,
    min_len: Option<usize>,
    total_min_len: Arc<Mutex<usize>>,

    phred: u8,
    max_qva: u8,
    total_max_qva: Arc<Mutex<u8>>,
    each: HashMap<usize, Vec<usize>>,
    total_each: Arc<Mutex<HashMap<usize, Vec<usize>>>>,
}

impl Info {
    fn new(phred: u8, max_qva: u8) -> Self {
        Info {
            num_a: 0,
            total_num_a: Arc::new(Mutex::new(0)),
            rate_a: 0.0,
            num_t: 0,
            total_num_t: Arc::new(Mutex::new(0)),
            rate_t: 0.0,
            num_g: 0,
            total_num_g: Arc::new(Mutex::new(0)),
            rate_g: 0.0,
            num_c: 0,
            total_num_c: Arc::new(Mutex::new(0)),
            rate_c: 0.0,
            num_n: 0,
            total_num_n: Arc::new(Mutex::new(0)),
            rate_n: 0.0,
            num_read: 0,
            total_num_read: Arc::new(Mutex::new(0)),
            num_base: 0,
            total_num_base: Arc::new(Mutex::new(0)),
            num_q5: 0,
            total_num_q5: Arc::new(Mutex::new(0)),
            rate_q5: 0.0,
            num_q10: 0,
            total_num_q10: Arc::new(Mutex::new(0)),
            rate_q10: 0.0,
            num_q15: 0,
            total_num_q15: Arc::new(Mutex::new(0)),
            rate_q15: 0.0,
            num_q20: 0,
            total_num_q20: Arc::new(Mutex::new(0)),
            rate_q20: 0.0,
            num_q30: 0,
            total_num_q30: Arc::new(Mutex::new(0)),
            rate_q30: 0.0,
            rate_gc: 0.0,
            ave_len: 0.0,
            max_len: 0,
            total_max_len: Arc::new(Mutex::new(0)),
            min_len: None,
            total_min_len: Arc::new(Mutex::new(0)),
            phred,
            max_qva,
            total_max_qva: Arc::new(Mutex::new(0)),
            each: HashMap::new(),
            total_each: Arc::new(Mutex::new(HashMap::new())),
        }
    }
    fn calc(&mut self) {
        let num_base = *self.total_num_base.lock();
        let num_read = *self.total_num_read.lock();
        self.ave_len = num_base as f64 / num_read as f64;
        self.rate_a = *self.total_num_a.lock() as f64 / num_base as f64;
        self.rate_t = *self.total_num_t.lock() as f64 / num_base as f64;
        self.rate_g = *self.total_num_g.lock() as f64 / num_base as f64;
        self.rate_c = *self.total_num_c.lock() as f64 / num_base as f64;
        self.rate_n = *self.total_num_n.lock() as f64 / num_base as f64;
        self.rate_gc =
            (*self.total_num_g.lock() + *self.total_num_c.lock()) as f64 / num_base as f64;
        self.rate_q5 = *self.total_num_q5.lock() as f64 / num_base as f64;
        self.rate_q10 = *self.total_num_q10.lock() as f64 / num_base as f64;
        self.rate_q15 = *self.total_num_q15.lock() as f64 / num_base as f64;
        self.rate_q20 = *self.total_num_q20.lock() as f64 / num_base as f64;
        self.rate_q30 = *self.total_num_q30.lock() as f64 / num_base as f64;
    }
}

impl ParallelProcessor for Info {
    fn process_record<Rf: Record>(&mut self, record: Rf) -> Result<(), ProcessError> {
        self.num_read += 1;
        self.num_base += record.seq().len();

        let this_q = *record.qual().unwrap().iter().max().unwrap() - self.phred;
        if this_q > self.max_qva {
            self.max_qva = this_q;
        }
        let len = record.seq().len();
        if len > self.max_len {
            self.max_len = len;
        }
        match self.min_len {
            Some(v) => {
                if v >= len {
                    self.min_len = Some(len)
                }
            }
            None => self.min_len = Some(len),
        }

        for &base in record.seq() {
            match base {
                b'A' => self.num_a += 1,
                b'T' => self.num_t += 1,
                b'G' => self.num_g += 1,
                b'C' => self.num_c += 1,
                b'N' => self.num_n += 1,
                _ => {}
            }
        }

        for (pos, (sf, sq)) in record
            .seq()
            .iter()
            .zip(record.qual_str().as_bytes())
            .enumerate()
        {
            let idx = (*sq - self.phred) as usize;
            if idx >= 5 {
                self.num_q5 += 1;
            }
            if idx >= 10 {
                self.num_q10 += 1;
            }
            if idx >= 15 {
                self.num_q15 += 1;
            }
            if idx >= 20 {
                self.num_q20 += 1;
            }
            if idx >= 30 {
                self.num_q30 += 1;
            }
            // for hashmap part
            if let std::collections::hash_map::Entry::Vacant(e) = self.each.entry(pos) {
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
                let tf = self.each.get_mut(&pos).unwrap();
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

        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<(), ProcessError> {

        *self.total_num_a.lock() += self.num_a;
        *self.total_num_t.lock() += self.num_t;
        *self.total_num_g.lock() += self.num_g;
        *self.total_num_c.lock() += self.num_c;
        *self.total_num_n.lock() += self.num_n;
        *self.total_num_read.lock() += self.num_read;
        *self.total_num_base.lock() += self.num_base;
        *self.total_num_q5.lock() += self.num_q5;
        *self.total_num_q10.lock() += self.num_q10;
        *self.total_num_q15.lock() += self.num_q15;
        *self.total_num_q20.lock() += self.num_q20;
        *self.total_num_q30.lock() += self.num_q30;

        *self.total_max_qva.lock() = self.max_qva;
        *self.total_max_len.lock() = self.max_len;
        *self.total_min_len.lock() = self.min_len.unwrap_or(0);

        let mut total_each = self.total_each.lock();
        for (pos, v) in self.each.iter() {
            if let std::collections::hash_map::Entry::Vacant(e) = total_each.entry(*pos) {
                e.insert(v.clone());
            } else {
                let tf = total_each.get_mut(pos).unwrap();
                for (i, val) in v.iter().enumerate() {
                    if i < tf.len() {
                        tf[i] += val;
                    } else {
                        tf.push(*val);
                    }
                }
            }
        }

        // reset for next batch
        self.num_a = 0;
        self.num_t = 0;
        self.num_g = 0;
        self.num_c = 0;
        self.num_n = 0;
        self.num_read = 0;
        self.num_base = 0;
        self.num_q5 = 0;
        self.num_q10 = 0;
        self.num_q15 = 0;
        self.num_q20 = 0;
        self.num_q30 = 0;
        self.max_len = 0;
        self.min_len = None;
        self.each = HashMap::new();
        Ok(())
    }
}

pub fn stat_fq(
    inp: Option<&String>,
    pre_sum: &String,
    pre_cyc: Option<&String>,
    phred: u8,
    ncp: usize,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    if ![33u8, 64u8].contains(&phred) {
        error!("{}", FqkitError::InvalidPhredValue);
        std::process::exit(1);
    }

    let fq_reader = fastq::Reader::new(file_reader(inp)?);
    info!("summary result write to file: {}", pre_sum);
    if let Some(file) = pre_cyc {
        info!("cycle result write to file: {}", file);
    } else {
        info!("cycle result write to stdout");
    }

    let mut stat = Info::new(phred, 0);
    let mut fo = file_writer(Some(pre_sum), compression_level, stdout_type)?;
    let mut fc = file_writer(pre_cyc, compression_level, stdout_type)?;
    fq_reader.process_parallel(stat.clone(), ncp)?;

    stat.calc();

    // output summary result
    writeln!(&mut fo, "read average length:\t{:.0}", stat.ave_len)?;
    writeln!(&mut fo, "read min length:\t{}", stat.total_min_len.lock())?;
    writeln!(&mut fo, "read max length:\t{}", stat.total_max_len.lock())?;
    writeln!(&mut fo, "total gc content(%):\t{:.2}", stat.rate_gc * 100.0)?;
    writeln!(&mut fo, "total read count:\t{}", stat.total_num_read.lock())?;
    writeln!(
        &mut fo,
        "total base count:\t{}\n",
        stat.total_num_base.lock()
    )?;
    writeln!(
        &mut fo,
        "base A count:\t{}\t({:.2}%)",
        stat.total_num_a.lock(),
        stat.rate_a * 100.0
    )?;
    writeln!(
        &mut fo,
        "base T count:\t{}\t({:.2}%)",
        stat.total_num_t.lock(),
        stat.rate_t * 100.0
    )?;
    writeln!(
        &mut fo,
        "base G count:\t{}\t({:.2}%)",
        stat.total_num_g.lock(),
        stat.rate_g * 100.0
    )?;
    writeln!(
        &mut fo,
        "base C count:\t{}\t({:.2}%)",
        stat.total_num_c.lock(),
        stat.rate_c * 100.0
    )?;
    writeln!(
        &mut fo,
        "base N count:\t{}\t({:.2}%)\n",
        stat.total_num_n.lock(),
        stat.rate_n * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 5 or higher (Q5+) (%)\t{}\t({:.2}%)",
        stat.total_num_q5.lock(),
        stat.rate_q5 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 10 or higher (Q10+) (%)\t{}\t({:.2}%)",
        stat.total_num_q10.lock(),
        stat.rate_q10 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 15 or higher (Q15+) (%)\t{}\t({:.2}%)",
        stat.total_num_q15.lock(),
        stat.rate_q15 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 20 or higher (Q20+) (%)\t{}\t({:.2}%)",
        stat.total_num_q20.lock(),
        stat.rate_q20 * 100.0
    )?;
    writeln!(
        &mut fo,
        "Number of base calls with quality value of 30 or higher (Q30+) (%)\t{}\t({:.2}%)",
        stat.total_num_q30.lock(),
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

    let max_qva = *stat.total_max_qva.lock();
    for i in 0..=max_qva {
        header.push(format!("{}", i));
    }
    fc.write_all(header.join("\t").as_bytes())?;
    fc.write_all(b"\n")?;

    let each = &stat.total_each.lock();
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
