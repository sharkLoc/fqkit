use clap::Parser;
use std::io::Result;

use colored::*;
use fastq::*;
use plot::*;

mod fastq;
mod plot;
mod utils;

#[derive(Parser, Debug)]
#[command(
    author = "size_t",
    version = "version 0.1.0",
    about = "fqkit: a simple program for fastq file manipulation",
    long_about = None
)]
struct Args {
    #[clap(subcommand)]
    command: Subcli,
}

#[derive(Parser, Debug)]
#[allow(non_camel_case_types)]
enum Subcli {
    /// subsample sequences from big fastq file.
    subfq {
        /// input fastq[.gz] file.
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
    /// summary for fastq format file
    stats {
        /// fastq input file name
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

    /// line plot for fastq quaility stats
    plot {
        /// input cycle result data
        #[arg(short = 'd', long = "data")]
        data: String,

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
        /// input fastq[.gz] file.
        input: Option<String>,

        /// output file name or write to stdout
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
}

fn main() -> Result<()> {
    let arg = Args::parse();
    match arg.command {
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
                        eprintln!(
                            "{}",
                            "[error]: opt -r used, fastq data can't from stdin.".red()
                        );
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
                        eprintln!(
                            "{}",
                            "[info]: need fastq file, if data from stdin, ignore this info.".red()
                        );
                        if out.is_some() {
                            select_fastq2(&None, num, seed, &Some(out.unwrap().as_str()))?;
                        } else {
                            select_fastq2(&None, num, seed, &None)?;
                        }
                    }
                }
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
            prefix,
            width,
            height,
            ylim,
            types,
        } => {
            let df = cycle_data(&Some(&data))?;
            let _x = plot_line(df, prefix, width, height, ylim, &types);
        }
        Subcli::stats {
            input,
            phred,
            sum,
            cyc,
        } => {
            if input.is_none() {
                println!("{}", "[info]: type --help for more information!".red());
                std::process::exit(1);
            }
            if cyc.is_some() {
                stat_fq(&Some(&input.unwrap()), &sum, &Some(&cyc.unwrap()), phred)?;
            } else {
                stat_fq(&Some(&input.unwrap()), &sum, &None, phred)?;
            }
        }
    }
    Ok(())
}
