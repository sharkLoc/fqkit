use crate::utils::*;
use anyhow::Result;
use bio::io::fastq;
use crossbeam::channel::bounded;
use log::*;

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
        Base {
            a: 0,
            t: 0,
            g: 0,
            c: 0,
            n: 0,
            read: 0,
        }
    }
}

pub fn size_fastq(
    fq: Option<&String>,
    ncpu: usize,
    chunk: usize,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<()> {

    let fq_reader = fastq::Reader::new(file_reader(fq)?);
    if let Some(inp) = fq {
        info!("reading from file: {}", inp);
    } else {
        info!("reading from stdin");
    }

    let mut chunk = chunk;
    if chunk == 0 {
        warn!(
            "read conut in chunk can't be: {}, changed to default value.",
            chunk
        );
        chunk = 5000;
    }

    if let Some(file) = out {
        info!("output result write to file: {}", file);
    } else {
        info!("output result write to stdout");
    }
    let mut fo = file_writer(out, compression_level, stdout_type)?;

    let mut base = Base::new();
    let mut bases = 0usize;

    if ncpu == 0 || ncpu == 1 {
        for rec in fq_reader.records().map_while(Result::ok) {
            base.read += 1;
            for nt in rec.seq().iter() {
                match *nt {
                    b'A' => base.a += 1,
                    b'T' => base.t += 1,
                    b'G' => base.g += 1,
                    b'C' => base.c += 1,
                    b'N' => base.n += 1,
                    _ => unreachable!(),
                }
            }
        }
        bases = base.a + base.t + base.g + base.c + base.n;
    } else {
        let (tx, rx) = bounded(5000);//unbounded();
        let mut fqiter = fq_reader.records();
        loop {
            let chunk: Vec<_> = fqiter.by_ref().take(chunk).map_while(Result::ok).collect();
            if chunk.is_empty() {
                break;
            }
            tx.send(chunk).unwrap();
        }
        drop(tx);

        crossbeam::scope(|s| {
            let (tx2, rx2) = bounded(5000); //unbounded();
            let _handles: Vec<_> = (0..ncpu)
                .map(|_| {
                    let rx_tmp = rx.clone();
                    let tx_tmp = tx2.clone();
                    s.spawn(move |_| {
                        for vrec in rx_tmp {
                            let mut base = Base::new();
                            for rec in vrec {
                                base.read += 1;
                                for nt in rec.seq().iter() {
                                    match *nt {
                                        b'A' => base.a += 1,
                                        b'T' => base.t += 1,
                                        b'G' => base.g += 1,
                                        b'C' => base.c += 1,
                                        b'N' => base.n += 1,
                                        _ => unreachable!(),
                                    }
                                }
                            }
                            tx_tmp.send(base).unwrap();
                        }
                    });
                })
                .collect();
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
        })
        .unwrap();
    }

    fo.write_all(
        format!(
            "reads:{}\tbases:{}\tA:{}\tT:{}\tG:{}\tC:{}\tN:{}\n",
            base.read, bases, base.a, base.t, base.g, base.c, base.n
        )
        .as_bytes(),
    )?;
    fo.flush()?;

    Ok(())
}
