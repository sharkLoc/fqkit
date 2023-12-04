use anyhow::{Error,Ok};
use clap::Parser;
use chrono::Local;
use env_logger::{Builder,fmt::Color};
use log::{error, warn, LevelFilter,Level};
use size::size_fastq;
use view::view_fq;
use std::io::Write;


mod view;
mod size;
mod reverse;
use reverse::*;
mod trimfq;
use trimfq::*;
mod gcplot;
use gcplot::*;
mod split2;
use split2::*;
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
mod flatten;
use flatten::*;
mod utils;


#[derive(Parser, Debug)]
#[command(
    author = "sharkLoc",
    version = "0.2.16",
    about = "A simple program for fastq file manipulation",
    long_about = None,
    next_line_help = false,
)]
#[command(help_template = 
    "{name}: {about}\n\nVersion: {version}\
    \nAuthors: {author} <mmtinfo@163.com>\
    \n\n{usage-heading} {usage}\n\n{all-args}\n"
)]
struct Args {
    #[clap(subcommand)]
    command: Subcli,
    /// be quiet and do not show extra information
    #[arg(short = 'q', long = "quiet", global= true, help_heading = Some("Global FLAGS"))]
    pub quiet: bool,
}

#[derive(Parser, Debug)]
#[allow(non_camel_case_types)]
enum Subcli {
    /// get first N records from fastq file
    topn {
        /// input fasta[.gz] file, or read from stdin
        input: Option<String>,
        /// print first N fasta records
        #[arg(short = 'n', long = "num", default_value_t = 10)]
        num: usize,
        /// output fasta[.gz] file name or write to stdout, files ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    }, 
    /// subsample sequences from big fastq file.
    #[command(visible_alias = "sample")]
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
        #[arg(short = 'r', long = "rdc", help_heading = Some("FLAGS"))]
        rdc: bool,
        /// fastq output file name or write to stdout, files ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// trim fastq file
    trim {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// trim int bp from left  
        #[arg(short, long, default_value_t=0)]
        left: usize,
        /// trim int bp from right
        #[arg(short, long, default_value_t=0)]
        right: usize,
        /// fastq output file name or write to stdout, files ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// search reads/motifs from fastq file
    search {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// specify uppercase pattern/motif, regular expression supported, e.g., -p "ATC{2,}" or -p "ATCCG"
        /// for multiple motifs, -p "TTAGGG|CCCTAA"
        #[arg(short = 'p', long = "pattern",verbatim_doc_comment)]
        pat: String,
        /// number of additional worker threads to use
        #[arg(short='@', long="thread", default_value_t = 1)]
        thread: usize,
        /// output contain pattern/motif reads result fastq[.gz] file or write to stdout,
        /// file ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out", verbatim_doc_comment)]
        out: Option<String>,
    },
    /// summary for fastq format file
    #[command(visible_alias = "stat", subcommand_help_heading = Some("Statistics"))]
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
    /// report the number sequences and bases
    size {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// number of additional worker threads to use
        #[arg(short='@', long="thread", default_value_t = 1)]
        thread: usize,
        /// output file name or write to stdout, file ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,

    },
    /// line plot for A T G C N percentage in read position
    plot {
        /// input cycle result data: fqkit stats cycle output
        #[arg(short = 'd', long = "data")]
        data: String,
        /// if specified, show line plot in terminal
        #[arg( short = 's', long ="show-terminal", help_heading = Some("FLAGS"))]
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
        /// if specified, remove sequence id description
        #[arg(short='r', long="remove", help_heading = Some("FLAGS"))]
        remove: bool,
        /// output file name or write to stdout, file ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// flatten fastq sequences
    flatten {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// filed number, id:1, sequence:2, symbol:4, quality:8
        /// eg. output id, sequence and quality value: 1 + 2 + 8 == 11 ,
        #[arg(short = 'f', long = "field", default_value_t = 3, verbatim_doc_comment)]
        flag: u8,
        /// output seprater, can be ",",  ";", 
        #[arg(short = 's', long = "sep", default_value_t='\t')]
        sep: char,
        /// output file name[.gz] or write to stdout, file ending in .gz will be compressed automatically
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
        #[arg(short = 'r', long = "rev_comp", help_heading = Some("FLAGS"))]
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
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// output file name[.gz] or write to stdout, file ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
        /// read name list file, one name per line and without read name prefix "@"
        #[arg(short = 'n', long = "name")]
        name: String,
        /// save removed reads in read name list
        #[arg(short = 's', long = "save",default_value_t=String::from("rm.fq.gz"))]
        save: String,
    },
    /// get a reverse-complement of fastq file.
    #[command(visible_alias = "rev")]
    reverse {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// if set, just output reverse sequences, the quality scores are also reversed
        #[arg(short = 'r', long = "reverse", help_heading = Some("FLAGS"))]
        rev: bool,
        /// output file name[.gz] or write to stdout, file ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// split interleaved fastq file
    split {
        /// input fastq[.gz] file, or read from stdin
        #[arg(short = 'i', long = "input")]
        input: Option<String>,
        /// output fastq file prefix name
        #[arg(short = 'p', long = "prefix")]
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
    /// split fastq file by records number
    split2 {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// set record number for each mini fastq file
        #[arg(short = 'n', long = "num", default_value_t = 200000)]
        num: usize,
        /// if specified, output gzip compressed file
        #[arg(short = 'z', long = "gzip", help_heading = Some("FLAGS"))]
        gzip: bool,
        /// output prefix name
        #[arg(short = 'p', long = "prefix", default_value_t = String::from("sub"))]
        name: String,
    },
    /// get GC content result and plot
    gcplot {
        /// input fastq[.gz] file, or read from stdin
        input: Option<String>,
        /// output GC contnet result file name
        #[arg(short = 'o', long = "out")]
        output: Option<String>,
        /// if specified, show histogram graphs in terminal
        #[arg( short = 's', long ="show-terminal", help_heading = Some("FLAGS"))]
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
    },
    /// view fastq file page by page
    view {
        /// input fastq[.gz] file
        input: Option<String>,
        /// output reads page by page, file name ending in .gz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
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
            if let Some(input) = input {
                if let Some(out) = out {
                    top_n_records(&Some(&input), num, &Some(&out), arg.quiet)?;
                } else {
                    top_n_records(&Some(&input), num, &None, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    top_n_records(&None, num, &Some(&out), arg.quiet)?;
                } else {
                    top_n_records(&None, num, &None, arg.quiet)?;
                }
            }
        }
        Subcli::subfq { input, seed, num, rdc, out,} => {
            if rdc {
                match input {
                    Some(x) => {
                        if out.is_some() {
                            select_fastq(&Some(x.as_str()),num, seed, &Some(out.unwrap().as_str()), arg.quiet )?;
                        } else {
                            select_fastq(&Some(x.as_str()), num, seed, &None, arg.quiet)?;
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
                            select_fastq2(&Some(x.as_str()),num, seed,&Some(out.unwrap().as_str()), arg.quiet)?;
                        } else {
                            select_fastq2(&Some(x.as_str()), num, seed, &None, arg.quiet)?;
                        }
                    }
                    None => {
                        warn!("need fastq file, if data from stdin, ignore this info.");
                        if out.is_some() {
                            select_fastq2(&None, num, seed, &Some(out.unwrap().as_str()), arg.quiet)?;
                        } else {
                            select_fastq2(&None, num, seed, &None, arg.quiet)?;
                        }
                    }
                }
            }
        }
        Subcli::trim { input, left, right, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    trim_fq(&Some(input.as_str()), left, right, &Some(out.as_str()), arg.quiet)?;
                } else {
                    trim_fq(&Some(input.as_str()), left, right, &None, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    trim_fq(&None, left, right, &Some(out.as_str()), arg.quiet)?;
                } else {
                    trim_fq(&None, left, right, &None, arg.quiet)?;
                }
            }
        }
        Subcli::search { input, pat, thread, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    search_fq(&Some(&input), &pat, &Some(&out), thread, arg.quiet)?;
                }else {
                    search_fq(&Some(&input), &pat, &None, thread, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    search_fq(&None, &pat, &Some(&out), thread, arg.quiet)?;
                }else {
                    search_fq(&None, &pat, &None, thread, arg.quiet)?;
                }
            }
        }
        Subcli::fq2fa { input, remove, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    fq2fa(&Some(&input), remove, &Some(&out), arg.quiet)?;
                } else {
                    fq2fa(&Some(&input), remove, &None, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    fq2fa(&None, remove, &Some(&out), arg.quiet)?;
                } else {
                    fq2fa(&None, remove, &None, arg.quiet)?;
                }
            }
        }
        Subcli::flatten { input, flag, sep, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    flatten_fq(&Some(&input), &Some(&out), flag, sep, arg.quiet)?;
                } else {
                    flatten_fq(&Some(&input), &None, flag, sep, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    flatten_fq(&None, &Some(&out), flag, sep, arg.quiet)?;
                } else {
                    flatten_fq(&None, &None, flag, sep, arg.quiet)?;
                }
            }
        },
        Subcli::plot { data, show, prefix, width, height, ylim,types,} => {
            let df = cycle_data(&Some(&data))?;
            let _x = plot_line(df, show, prefix, width, height, ylim, &types, arg.quiet);
        }
        Subcli::stats { input, phred, sum,cyc,} => {
            if let Some(input) = input {
                if let Some(cyc) = cyc {
                    stat_fq(&Some(&input), &sum, &Some(&cyc), phred, arg.quiet)?;
                } else {
                    stat_fq(&Some(&input), &sum, &None, phred, arg.quiet)?;
                }
            } else {
                if let Some(cyc) = cyc {
                    stat_fq(&None, &sum, &Some(&cyc), phred, arg.quiet)?;
                } else {
                    stat_fq(&None, &sum, &None, phred, arg.quiet)?;
                }
            }
        }
        Subcli::size { input, thread, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    size_fastq(&Some(&input), thread, &Some(&out), arg.quiet)?;
                } else {
                    size_fastq(&Some(&input), thread, &None, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    size_fastq(&None, thread, &Some(&out), arg.quiet)?;
                } else {
                    size_fastq(&None, thread, &None, arg.quiet)?;
                }
            }

        }
        Subcli::barcode { read1, read2, bar, mode, trans, mismatch, outdir, } => {
               split_fq(&read1, &read2, &bar, trans, mode, mismatch, &outdir, arg.quiet)?; 
        }
        Subcli::remove { input, out, name , save} => {
            if let Some(input) = input {
                if let Some(out) = out {
                    remove_read(&Some(&input), &Some(&out) ,&name, &save, arg.quiet)?;
                } else {
                    remove_read(&Some(&input), &None ,&name, &save, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    remove_read(&None, &Some(&out) ,&name, &save, arg.quiet)?;
                } else {
                    remove_read(&None, &None ,&name, &save, arg.quiet)?;
                }
            }
        }
        Subcli::reverse { input, rev, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    reverse_comp_seq(&Some(&input), &Some(&out), rev, arg.quiet)?;
                } else {
                    reverse_comp_seq(&Some(&input), &None, rev, arg.quiet)?;
                }   
            } else {
                if let Some(out) = out {
                    reverse_comp_seq(&None, &Some(&out), rev, arg.quiet)?;
                } else {
                    reverse_comp_seq(&None, &None, rev, arg.quiet)?;
                }
            }
        }
        Subcli::split { input, pre,   out, } => {
            if let Some(input) = input {
                split_interleaved(&Some(&input), &out, &pre, arg.quiet)?; 
            } else {
                split_interleaved(&None, &out, &pre, arg.quiet)?; 
            }
        }
        Subcli::merge { read1,   read2,  out, } => {
            interleaved(&Some(read1.as_str()), &Some(read2.as_str()), &Some(out.as_str()), arg.quiet)?;    
        }
        Subcli::split2 { input, num, gzip, name } => {
            if let Some(input) = input {
                split_chunk(&Some(&input), num, gzip,&name, arg.quiet)?;
            } else {
                split_chunk(&None, num, gzip,&name, arg.quiet)?;
            }
        }
        Subcli::gcplot { input, output, show, prefix, width, height, ylim, types } => {
            if let Some(input) = input {
                if let Some(output) = output {
                    gc_content(&Some(&input), &Some(&output), show, prefix, width, height, ylim, &types, arg.quiet)?;
                } else {
                    gc_content(&Some(&input),&None, show, prefix, width, height, ylim, &types, arg.quiet)?;
                }
            } else {
                if let Some(output) = output {
                    gc_content(&None,&Some(&output), show, prefix, width, height, ylim, &types, arg.quiet)?;
                } else {
                    gc_content(&None,&None, show, prefix, width, height, ylim, &types, arg.quiet)?;
                }
            }
        }
        Subcli::view { input, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    view_fq(&Some(&input), &Some(&out), arg.quiet)?;
                } else {
                    view_fq(&Some(&input), &None, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    view_fq(&None, &Some(&out), arg.quiet)?;
                } else {
                    view_fq(&None, &None, arg.quiet)?;
                }
            }
        }
    }

    Ok(())
}
