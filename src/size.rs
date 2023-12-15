use bio::io::fastq;
use std:: io::Result;
use log::*;
use crossbeam::channel::unbounded;
use crate::utils::*;
use std::time::Instant;

#[derive(Clone, Copy)]
struct Base {
    a: usize,
    t: usize,
    g: usize,
    c: usize,
    n: usize,
    read: usize,
}

impl Base {
    pub fn new() -> Self {
        Base { a: 0, t: 0, g: 0, c: 0, n: 0, read: 0 }
    }
}

pub fn size_fastq(
    fq: &Option<&str>,
    ncpu: usize,
    chunk: usize,
    out: &Option<&str>, 
    quiet: bool,
) -> Result<()> {
    let start = Instant::now();
    if !quiet {
        if let Some(inp) = fq {
            info!("reading from file: {}", inp);
        } else {
            info!("reading from stdin");
        }
        info!("additional worker threads is: {}", ncpu);
        if let Some(file) = out {
            info!("output result write to file: {}", file);
        } else {
            info!("output result write to stdout");
        }
    }

    let mut chunk = chunk;
    if chunk == 0 {
        warn!("read conut in chunk can't be: {}, changed to default value.",chunk);
        chunk = 5000;
    }
    let fq = fastq::Reader::new(file_reader(fq)?);
    let mut fo = file_writer(out)?;
    let mut base = Base::new();
    let mut bases = 0usize;

    if ncpu == 0 || ncpu == 1 {
        for rec in fq.records().flatten() {
            base.read += 1;
            for nt in rec.seq().iter() {
                match nt {
                    &b'A' => base.a +=1,
                    &b'T' => base.t +=1,
                    &b'G' => base.g +=1,
                    &b'C' => base.c +=1,
                    &b'N' => base.n +=1,
                    _ => unreachable!(),
                }
            };
        }
        bases = base.a + base.t + base.g + base.c + base.n; 
    } else {
        let (tx, rx) = unbounded();
        let mut fqiter = fq.records();
        loop {
            let chunk: Vec<_> = fqiter.by_ref().take(chunk).flatten().collect();
            if chunk.is_empty() { break; }
            tx.send(chunk).unwrap();
        }
        drop(tx);
    
        crossbeam::scope(|s| {
            let (tx2, rx2) = unbounded();
            let _handles: Vec<_> = (0..ncpu).map(|_| {
                let rx_tmp = rx.clone();
                let tx_tmp = tx2.clone();
                s.spawn(move |_| {
                    for vrec in rx_tmp {
                        let mut base = Base::new();
                        for rec in vrec {
                            base.read += 1;
                            for nt in rec.seq().iter() {
                                match nt {
                                    &b'A' => base.a +=1,
                                    &b'T' => base.t +=1,
                                    &b'G' => base.g +=1,
                                    &b'C' => base.c +=1,
                                    &b'N' => base.n +=1,
                                    _ => unreachable!(),
                                }
                            }
                        };
                        tx_tmp.send(base).unwrap();
                    }
                });
          
            }).collect();
            drop(tx2);
    
            for data in rx2.iter() {
                base.read += data.read;
                base.a += data.a;
                base.t += data.t;
                base.g += data.g;
                base.c += data.c;
                base.n += data.n;
            }
            bases = base.a + base.t + base.g + base.c + base.n; 
        }).unwrap();      
    }

    fo.write(format!("reads:{}\tbases:{}\tA:{}\tT:{}\tG:{}\tC:{}\tN:{}\n", base.read,bases,base.a,base.t,base.g,base.c,base.n).as_bytes()).unwrap();
    fo.flush()?;
    
    if !quiet {
        info!("time elapsed is: {:?}",start.elapsed());
    }
    Ok(())
}
