use anyhow::{Error,Ok};
use clap::Parser;
use chrono::Local;
use env_logger::{Builder,fmt::Color};
use log::{error, warn, info, LevelFilter,Level};
use std::io::Write;

mod gcplot;
use gcplot::*;
mod split;
use split::*;
mod search;
use search::*;
mod subfq;
use subfq::*;
mod merge;
use merge::*;
mod remove;
use remove::*;
mod stats;
use stats::*;
mod top;
use top::*;
mod plot;
use plot::*;
mod barcode;
use barcode::*;
mod fq2fa;
use fq2fa::*;
mod utils;


#[derive(Parser, Debug)]
#[command(
    author = "size_t",
    version = "version 0.2.8",
    about = "fqkit: a simple program for fastq file manipulation",
    long_about = None,
    next_line_help = true
)]
struct Args {
    #[clap(subcommand)]
    command: Subcli,
}

#[derive(Parser, Debug)]
#[allow(non_camel_case_types)]
enum Subcli {
    /// get first N records from fastq file
    topn {
        /// input fasta[.gz] file
        #[arg(short = 'i', long = "input")]
        input: String,
        /// print first N fasta records
        #[arg(short = 'n', long = "num", default_value_t = 10)]
        num: usize,
        /// output fasta[.gz] file name, or write to stdout
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    }, 
    /// subsample sequences from big fastq file.
    subfq {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// set rand seed.
        #[arg(short = 's', long = "seed", default_value_t = 69)]
        seed: u64,
        /// subseq number
        #[arg(short = 'n', long = "num")]
        num: usize,
        /// reduce much memory but cost more time
        #[arg(short = 'r', long = "rdc")]
        rdc: bool,
        /// read output file name or write to stdout
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// search reads/motifs from fastq file
    search {
        /// input fastq[.gz] file, or read from stdin
        input: String,
        /// specify uppercase pattern/motif, e.g., -p "ATC{2,}" or -p ATCCG
        #[arg(short = 'p', long = "pattern")]
        pat: String,
        /// output contain pattern/motif reads result fastq[.gz] file, or write to stdout
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// summary for fastq format file
    stats {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        ///phred score 33 or 64
        #[arg(short = 'p', long = "phred", default_value_t = 33)]
        phred: u8,
        /// specify a name for summary output file
        #[arg(short='s',long="sumy",default_value_t=String::from("summary.txt"))]
        sum: String,
        /// if not specified, cycle result write to stdout
        #[arg(short = 'c', long = "cycle")]
        cyc: Option<String>,
    },
    /// line plot for A T G C N percentage in read position
    plot {
        /// input cycle result data: fqkit stats cycle output
        #[arg(short = 'd', long = "data")]
        data: String,
        /// if specified, show line plot in terminal
        #[arg( short = 's', long ="show-terminal")]
        show: bool,
        /// output base figure prefix name
        #[arg(short='p', long="prefix", default_value_t=String::from("base_plot"))]
        prefix: String,
        /// set output figure width
        #[arg(short = 'W', long = "width", default_value_t = 960)]
        width: usize,
        /// set output figure height
        #[arg(short = 'H', long = "height", default_value_t = 540)]
        height: usize,
        /// set max ylim (0~100)
        #[arg(short = 'y', long = "ylim", default_value_t = 50.0)]
        ylim: f32,
        /// figure type 'png' or 'svg'
        #[arg(short='t', long="types", default_value_t=String::from("png"))]
        types: String,
    },
    /// translate fastq to fasta
    fq2fa {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// output file name or write to stdout
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// split barcode for PE reads
    barcode {
        /// input read1 fastq[.gz] file
        #[arg(short = '1', long = "read1")]
        read1: String,
        /// input read2 fastq[.gz] file <barcode in this file>
        #[arg(short = '2', long = "read2")]
        read2: String,
        /// barcode list file, format eg:
        /// ATGCAGTG    sample1
        /// TGCAGTAC    sample2
        #[arg(short = 'b', long = "barcode", verbatim_doc_comment)] 
        bar: String,
        /// barcode position mode, 1:left, 2:right
        #[arg(short = 'm', long = "mode", default_value_t = 2)]
        mode: usize,
        /// barcode reverse complement
        #[arg(short = 'r', long = "rev_comp")]
        trans: bool,
        /// barcode mismatch base count
        #[arg(short = 'e', long = "error", default_value_t = 0)]
        mismatch: usize,
        /// fastq file output dir.
        #[arg(short = 'o', long = "outdir", default_value_t = String::from("."))]
        outdir: String,
    },
    /// remove reads by read name.
    remove {
        /// input fastq[.gz] file.
        #[arg(short = 'i', long = "input")]
        input: String,
        /// output file name[.gz] or write to stdout
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
        /// read name list file, one name per line and without read name prefix "@"
        #[arg(short = 'n', long = "name")]
        name: String,
    },
    /// split interleaved fastq file
    split {
        /// input fastq[.gz] file.
        #[arg(short = 'i', long = "input")]
        input: String,
        /// output fastq file prefix name
        #[arg(short = 'p', long = "pre")]
        pre: String,
        /// fastq file outdir
        #[arg(short = 'o', long = "out", default_value_t = String::from("."))]
        out: String,
    },
    /// merge PE reads as interleaved fastq file
    merge {
        /// input read1 fastq[.gz] file.
        #[arg(short = '1', long = "read1")]
        read1: String,
        /// input read2 fastq[.gz] file.
        #[arg(short = '2', long = "read2")]
        read2: String,
        /// output interleaved fastq file name.
        #[arg(short = 'o', long = "out", default_value_t = String::from("interleaved.fq.gz"))]
        out: String,
    },
    /// get GC content result and plot
    gcplot {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// output GC contnet result file name
        #[arg(short = 'o', long = "out")]
        output: Option<String>,
        /// if specified, show histogram graphs in terminal
        #[arg( short = 's', long ="show-terminal")]
        show: bool,
        /// output base figure prefix name
        #[arg(short='p', long="prefix", default_value_t=String::from("gc_plot"))]
        prefix: String,
        /// set output figure width
        #[arg(short = 'W', long = "width", default_value_t = 960)]
        width: usize,
        /// set output figure height
        #[arg(short = 'H', long = "height", default_value_t = 540)]
        height: usize,
        /// set max ylim (0~100)
        #[arg(short = 'y', long = "ylim", default_value_t = 15)]
        ylim: usize,
        /// figure type 'png' or 'svg'
        #[arg(short='t', long="types", default_value_t=String::from("png"))]
        types: String,
    }
}



fn main() -> Result<(), Error> {
    let mut builder = Builder::from_default_env();
    builder.format(|buf, record| {
        let mut style = buf.style();
        match record.level() {
            Level::Error => {
                style.set_color(Color::Red).set_bold(true);
            }
            Level::Warn => {
                style.set_color(Color::Yellow).set_bold(true);
            }
            Level::Info => {
                style.set_color(Color::Green).set_bold(true);
            }
            Level::Debug => {
                style.set_color(Color::Blue).set_bold(true);
            }
            Level::Trace => {
                style.set_color(Color::Magenta).set_bold(true);
            }
        }
        writeln!(buf,
            "[{} {} - {}] {}",
            Local::now().format("%Y-%m-%dT%H:%M:%S"),
            style.value(record.level()),
            buf.style().set_color(Color::Rgb(90, 150, 150)).value(record.target()),
            record.args()
        )
    })
    .filter(None, LevelFilter::Info)
    .init();
   

    let arg = Args::parse();
    match arg.command {
        Subcli::topn { input, num, out } => {
            if let Some(out) = out {
                top_n_records(&Some(&input), num, &Some(&out))?;
            } else {
                top_n_records(&Some(&input), num, &None)?;
            }
        }
        Subcli::subfq {
            input,
            seed,
            num,
            rdc,
            out,
        } => {
            if rdc {
                match input {
                    Some(x) => {
                        if out.is_some() {
                            select_fastq(
                                &Some(x.as_str()),
                                num,
                                seed,
                                &Some(out.unwrap().as_str()),
                            )?;
                        } else {
                            select_fastq(&Some(x.as_str()), num, seed, &None)?;
                        }
                    }
                    None => {
                        error!("opt -r used, fastq data can't from stdin.");
                        std::process::exit(1);
                    }
                }
            } else {
                match input {
                    Some(x) => {
                        if out.is_some() {
                            select_fastq2(
                                &Some(x.as_str()),
                                num,
                                seed,
                                &Some(out.unwrap().as_str()),
                            )?;
                        } else {
                            select_fastq2(&Some(x.as_str()), num, seed, &None)?;
                        }
                    }
                    None => {
                        warn!("need fastq file, if data from stdin, ignore this info.");
                        if out.is_some() {
                            select_fastq2(&None, num, seed, &Some(out.unwrap().as_str()))?;
                        } else {
                            select_fastq2(&None, num, seed, &None)?;
                        }
                    }
                }
            }
        }
        Subcli::search { input, pat, out } => {
            if let Some(out) = out {
                search_fq(&input, &pat, &Some(&out))?;
            }else {
                search_fq(&input, &pat, &None)?;
            }
        }
        Subcli::fq2fa { input, out } => match input {
            Some(x) => {
                if let Some(o) = out {
                    fq2fa(&Some(&x), &Some(&o))?;
                } else {
                    fq2fa(&Some(&x), &None)?;
                }
            }
            None => {
                if let Some(o) = out {
                    fq2fa(&None, &Some(&o))?;
                } else {
                    fq2fa(&None, &None)?;
                }
            }
        },
        Subcli::plot {
            data,
            show,
            prefix,
            width,
            height,
            ylim,
            types,
        } => {
            let df = cycle_data(&Some(&data))?;
            let _x = plot_line(df, show, prefix, width, height, ylim, &types);
        }
        Subcli::stats {
            input,
            phred,
            sum,
            cyc,
        } => {
            if input.is_none() {
                info!("type --help for more information!");
                std::process::exit(1);
            }
            if cyc.is_some() {
                stat_fq(&Some(&input.unwrap()), &sum, &Some(&cyc.unwrap()), phred)?;
            } else {
                stat_fq(&Some(&input.unwrap()), &sum, &None, phred)?;
            }
        }
        Subcli::barcode {
            read1,
            read2,
            bar,
            mode,
            trans,
            mismatch,
            outdir,
        } => {
               split_fq(&read1, &read2, &bar, trans, mode, mismatch, &outdir)?; 
        }
        Subcli::remove {
            input,
            out,
            name,
        } => {
            if out.is_some() {
                remove_read(&Some(&input), &Some(&out.unwrap()) ,&name)?;
            } else {
                remove_read(&Some(&input), &None ,&name)?;     
            }
        }
        Subcli::split {
            input,
            pre, 
            out,    
        } => {
            split_interleaved(&Some(&input), &out, &pre)?;    
        }
        Subcli::merge {
            read1, 
            read2,
            out,
        } => {
            interleaved(&Some(read1.as_str()), &Some(read2.as_str()), &Some(out.as_str()))?;    
        }
        Subcli::gcplot { input, output, show, prefix, width, height, ylim, types } => {
            if let Some(input) = input {
                if let Some(output) = output {
                    gc_content(&Some(&input), &Some(&output), show, prefix, width, height, ylim, &types)?;
                } else {
                    gc_content(&Some(&input),&None, show, prefix, width, height, ylim, &types)?;
                }
            } else {
                if let Some(output) = output {
                    gc_content(&None,&Some(&output), show, prefix, width, height, ylim, &types)?;
                } else {
                    gc_content(&None,&None, show, prefix, width, height, ylim, &types)?;
                }
            }
        }
    }

    Ok(())
}
