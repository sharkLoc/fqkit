use anyhow::{Error,Ok};
use clap::Parser;
use log::{error, warn};


mod slide;
use slide::*;
mod fqscore;
use fqscore::*;
mod grep;
use grep::*;
mod concat;
use concat::*;
mod shuffle;
use shuffle::*;
mod check;
use check::*;
mod mask;
use mask::*;
mod range;
use range::*;
mod sort;
use sort::*;
mod view;
use view::*;
mod size;
use size::*;
mod fq2sam;
use fq2sam::*;
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
mod command;
use command::*;
mod logger;
use logger::*;
mod utils;


fn main() -> Result<(), Error> {
    
    let arg = Args::parse();
    match arg.logfile {
        Some(v) => logger(arg.verbose, &Some(&v))? , 
        None => logger(arg.verbose, &None)?
    }
   
    match arg.command {
        Subcli::topn { input, num, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    top_n_records(&Some(&input), num, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    top_n_records(&Some(&input), num, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    top_n_records(&None, num, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    top_n_records(&None, num, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::subfq { input, seed, num, rdc, out,} => {
            if rdc {
                match input {
                    Some(x) => {
                        if out.is_some() {
                            select_fastq(&Some(x.as_str()),num, seed, &Some(out.unwrap().as_str()), arg.compression_level, arg.quiet )?;
                        } else {
                            select_fastq(&Some(x.as_str()), num, seed, &None, arg.compression_level, arg.quiet)?;
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
                            select_fastq2(&Some(x.as_str()),num, seed,&Some(out.unwrap().as_str()), arg.compression_level, arg.quiet)?;
                        } else {
                            select_fastq2(&Some(x.as_str()), num, seed, &None, arg.compression_level, arg.quiet)?;
                        }
                    }
                    None => {
                        warn!("need fastq file, if data from stdin, ignore this info.");
                        if out.is_some() {
                            select_fastq2(&None, num, seed, &Some(out.unwrap().as_str()), arg.compression_level, arg.quiet)?;
                        } else {
                            select_fastq2(&None, num, seed, &None, arg.compression_level, arg.quiet)?;
                        }
                    }
                }
            }
        }
        Subcli::trim { input, left, right, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    trim_fq(&Some(input.as_str()), left, right, &Some(out.as_str()), arg.compression_level, arg.quiet)?;
                } else {
                    trim_fq(&Some(input.as_str()), left, right, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    trim_fq(&None, left, right, &Some(out.as_str()), arg.compression_level, arg.quiet)?;
                } else {
                    trim_fq(&None, left, right, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::range { input, skip, take, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    range_fastq(&Some(&input), skip, take, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    range_fastq(&Some(&input), skip, take, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    range_fastq(&None, skip, take, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    range_fastq(&None, skip, take, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::search { input, pat, case, chunk, thread, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    search_fq(&Some(&input), &pat, case, chunk,&Some(&out), thread, arg.compression_level, arg.quiet)?;
                }else {
                    search_fq(&Some(&input), &pat, case, chunk,&None, thread, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    search_fq(&None, &pat, case, chunk,&Some(&out), thread, arg.compression_level, arg.quiet)?;
                }else {
                    search_fq(&None, &pat, case, chunk, &None, thread, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::grep { input, ids, full, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    grep_fastq(&Some(&input), &ids, full, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    grep_fastq(&Some(&input), &ids, full, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    grep_fastq(&None, &ids, full, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    grep_fastq(&None, &ids, full, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::fq2fa { input, remove, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    fq2fa(&Some(&input), remove, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    fq2fa(&Some(&input), remove, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    fq2fa(&None, remove, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    fq2fa(&None, remove, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::fq2sam { r1, r2, sm, rg, lb, pl, out } => {
            if let Some(r2) = r2 {
                if let Some(out) = out {
                    fastq2sam(&r1,&Some(&r2),&sm,rg,lb,pl,&Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    fastq2sam(&r1,&Some(&r2),&sm,rg,lb,pl,&None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    fastq2sam(&r1,&None,&sm,rg,lb,pl,&Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    fastq2sam(&r1,&None,&sm,rg,lb,pl,&None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::fqscore { input, to33, to64, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    phred_score(&Some(&input), &Some(&out), to33, to64, arg.compression_level, arg.quiet)?;
                } else {
                    phred_score(&Some(&input), &None, to33, to64, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    phred_score(&None, &Some(&out), to33, to64, arg.compression_level, arg.quiet)?;
                 } else {
                    phred_score(&None, &None, to33, to64, arg.compression_level, arg.quiet)?;
                 }
            }
        }
        Subcli::flatten { input, flag, sep, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    flatten_fq(&Some(&input), &Some(&out), flag, sep, arg.compression_level, arg.quiet)?;
                } else {
                    flatten_fq(&Some(&input), &None, flag, sep, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    flatten_fq(&None, &Some(&out), flag, sep, arg.compression_level, arg.quiet)?;
                } else {
                    flatten_fq(&None, &None, flag, sep, arg.compression_level, arg.quiet)?;
                }
            }
        },
        Subcli::plot { data, show, prefix, width, height, ylim,types,} => {
            let df = cycle_data(&Some(&data))?;
            let _x = plot_line(df, show, prefix, width, height, ylim, &types, arg.quiet);
        }
        Subcli::check { input, save, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    check_fastq(&Some(&input), save, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    check_fastq(&Some(&input), save, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    check_fastq(&None, save, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    check_fastq(&None, save, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::stats { input, phred, sum,cyc,} => {
            if let Some(input) = input {
                if let Some(cyc) = cyc {
                    stat_fq(&Some(&input), &sum, &Some(&cyc), phred, arg.compression_level, arg.quiet)?;
                } else {
                    stat_fq(&Some(&input), &sum, &None, phred, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(cyc) = cyc {
                    stat_fq(&None, &sum, &Some(&cyc), phred, arg.compression_level, arg.quiet)?;
                } else {
                    stat_fq(&None, &sum, &None, phred, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::shuffle { input, seed, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    shuffle_fastq(&Some(&input), seed, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    shuffle_fastq(&Some(&input), seed, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    shuffle_fastq(&None, seed, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    shuffle_fastq(&None, seed, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::size { input, thread, chunk, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    size_fastq(&Some(&input), thread, chunk, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    size_fastq(&Some(&input), thread, chunk, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    size_fastq(&None, thread, chunk, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    size_fastq(&None, thread, chunk, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::slide { input, window, step, suffix, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    slide_fastq(&Some(&input), step, window, &Some(&out), &suffix, arg.compression_level, arg.quiet)?;
                } else {
                    slide_fastq(&Some(&input), step, window, &None, &suffix, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    slide_fastq(&None, step, window, &Some(&out), &suffix, arg.compression_level, arg.quiet)?;
                } else {
                    slide_fastq(&None, step, window, &None, &suffix, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::sort { input, name, seq, gc, length, reverse, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    sort_fastq(&Some(&input), name, seq, gc, length, reverse, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    sort_fastq(&Some(&input), name, seq, gc, length, reverse, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    sort_fastq(&None, name, seq, gc, length, reverse, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    sort_fastq(&None, name, seq, gc, length, reverse, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::barcode { read1, read2, bar, mode, trans, mismatch, outdir, } => {
               split_fq(&read1, &read2, &bar, trans, mode, mismatch, &outdir, arg.compression_level, arg.quiet)?; 
        }
        Subcli::concat { read1, read2, out1, out2 } => {
            concat_fqstq_lane(&read1, &read2, &out1, &out2, arg.compression_level, arg.quiet)?;
        }
        Subcli::remove { input, out, name , save} => {
            if let Some(input) = input {
                if let Some(out) = out {
                    remove_read(&Some(&input), &Some(&out) ,&name, &save, arg.compression_level, arg.quiet)?;
                } else {
                    remove_read(&Some(&input), &None ,&name, &save, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    remove_read(&None, &Some(&out) ,&name, &save, arg.compression_level, arg.quiet)?;
                } else {
                    remove_read(&None, &None ,&name, &save, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::reverse { input, rev, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    reverse_comp_seq(&Some(&input), &Some(&out), rev, arg.compression_level, arg.quiet)?;
                } else {
                    reverse_comp_seq(&Some(&input), &None, rev, arg.compression_level, arg.quiet)?;
                }   
            } else {
                if let Some(out) = out {
                    reverse_comp_seq(&None, &Some(&out), rev, arg.compression_level, arg.quiet)?;
                } else {
                    reverse_comp_seq(&None, &None, rev, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::split { input, pre,   out, } => {
            if let Some(input) = input {
                split_interleaved(&Some(&input), &out, &pre, arg.compression_level, arg.quiet)?; 
            } else {
                split_interleaved(&None, &out, &pre, arg.compression_level, arg.quiet)?; 
            }
        }
        Subcli::merge { read1,   read2,  out, } => {
            interleaved(&Some(read1.as_str()), &Some(read2.as_str()), &Some(out.as_str()), arg.compression_level, arg.quiet)?;    
        }
        Subcli::mask { input, phred, low, chars, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    mask_fastq(&Some(&input), phred, low, chars, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    mask_fastq(&Some(&input), phred, low, chars, &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    mask_fastq(&None, phred, low, chars, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    mask_fastq(&None, phred, low, chars, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::split2 { input, num, gzip, name } => {
            if let Some(input) = input {
                split_chunk(&Some(&input), num, gzip,&name, arg.compression_level, arg.quiet)?;
            } else {
                split_chunk(&None, num, gzip,&name, arg.compression_level, arg.quiet)?;
            }
        }
        Subcli::gcplot { input, output, show, prefix, width, height, ylim, types } => {
            if let Some(input) = input {
                if let Some(output) = output {
                    gc_content(&Some(&input), &Some(&output), show, prefix, width, height, ylim, &types, arg.compression_level, arg.quiet)?;
                } else {
                    gc_content(&Some(&input),&None, show, prefix, width, height, ylim, &types, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(output) = output {
                    gc_content(&None,&Some(&output), show, prefix, width, height, ylim, &types, arg.compression_level, arg.quiet)?;
                } else {
                    gc_content(&None,&None, show, prefix, width, height, ylim, &types, arg.compression_level, arg.quiet)?;
                }
            }
        }
        Subcli::view { input, out } => {
            if let Some(input) = input {
                if let Some(out) = out {
                    view_fq(&Some(&input), &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    view_fq(&Some(&input), &None, arg.compression_level, arg.quiet)?;
                }
            } else {
                if let Some(out) = out {
                    view_fq(&None, &Some(&out), arg.compression_level, arg.quiet)?;
                } else {
                    view_fq(&None, &None, arg.compression_level, arg.quiet)?;
                }
            }
        }
    }

    Ok(())
}