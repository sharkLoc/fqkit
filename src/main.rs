use anyhow::{Error,Ok};
use clap::Parser;


mod logger;
use logger::*;
mod utils;
mod command;
use command::*;
mod cli;
use cli::{length::*, rename::*, select::*, tail::*, top::*, filter::*, slide::*, flatten::*, fq2fa::*, fq2sam::*, fqscore::*,
    grep::*, concat::*, shuffle::*, check::*, mask::*, range::*, sort::*, view::*, size::*, reverse::*, trimfq::*, gcplot::*, 
    split::*, split2::*, search::*, subfq::*, merge::*, remove::*, stats::*, plot::*, barcode::*, 
};

fn main() -> Result<(), Error> {

    let arg = Args::parse();
    match arg.logfile {
        Some(v) => logger(arg.verbose, &Some(&v), arg.quiet)?,
        None => logger(arg.verbose, &None, arg.quiet)?
    }
   
    match arg.command {
        Subcli::topn { input, num, out } => {
            top_n_records(input.as_ref(), num, out.as_ref(), arg.compression_level)?;
        }
        Subcli::tail { input, num, out } => {
            tail_n_records(input.as_ref(), num, out.as_ref(), arg.compression_level)?;
        }
        Subcli::subfq { input, seed, num, rdc, out,} => {
            subset_fastq(rdc, input.as_ref(), num, seed, out.as_ref(), arg.compression_level)?;
        }
        Subcli::select { read1, read2, out1, out2 } => {
            select_pe_fastq(&read1, &read2, &out1, &out2, arg.compression_level)?;
        }
        Subcli::trim { input, left, right, len, out } => {
            trim_fq(input.as_ref(), left, right, len, out.as_ref(), arg.compression_level)?;
        }
        Subcli::range { input, skip, take, out } => {
            range_fastq(input.as_ref(), skip, take, out.as_ref(), arg.compression_level)?;
        }
        Subcli::search { input, pat, case, chunk, thread, out } => {
            search_fq(input.as_ref(), &pat, case, chunk, out.as_ref(), thread, arg.compression_level)?;
        }
        Subcli::grep { input, ids, full, out } => {
            grep_fastq(input.as_ref(), &ids, full, out.as_ref(), arg.compression_level)?;
        }
        Subcli::fq2fa { input, remove, out } => {
            fq2fa(input.as_ref(), remove, out.as_ref(), arg.compression_level)?;
        }
        Subcli::fq2sam { r1, r2, sm, rg, lb, pl, out } => {
            fastq2sam(&r1, r2.as_ref(), &sm,  rg, lb, pl, out.as_ref(), arg.compression_level)?;
        }
        Subcli::fqscore { input, to33, to64, out } => {
            phred_score(input.as_ref(), out.as_ref(), to33, to64, arg.compression_level)?;
        }
        Subcli::flatten { input, flag, sep, out } => {
            flatten_fq(input.as_ref(), out.as_ref(), flag, sep, arg.compression_level)?;
        },
        Subcli::plot { data, show, prefix, width, height, ylim,types,} => {
            let df = cycle_data(Some(&data))?;
            let _x = plot_line(df, show, prefix, width, height, ylim, &types);
        }
        Subcli::check { input, save, out } => {
            check_fastq(input.as_ref(), save, out.as_ref()  , arg.compression_level)?;
        }
        Subcli::stats { input, phred, sum,cyc,} => {
            stat_fq(input.as_ref(), &sum, cyc.as_ref(), phred, arg.compression_level)?;
        }
        Subcli::shuffle { input, seed, out } => {
            shuffle_fastq(input.as_ref(), seed, out.as_ref(), arg.compression_level)?;
        }
        Subcli::size { input, thread, chunk, out } => {
            size_fastq(input.as_ref(), thread, chunk, out.as_ref(), arg.compression_level)?;
        }
        Subcli::slide { input, window, step, suffix, out } => {
            slide_fastq(input.as_ref(), step, window, out.as_ref(), &suffix, arg.compression_level)?;
        }
        Subcli::sort { input, name, seq, gc, length, reverse, out } => {
            sort_fastq(input.as_ref(), name, seq, gc, length, reverse, out.as_ref(), arg.compression_level)?;
        }
        Subcli::barcode { read1, read2, bar, mode, trans, mismatch, gzip, bzip2, xz, outdir, } => {
               split_fq(&read1, &read2, &bar, trans, mode, mismatch, &outdir, gzip, bzip2, xz,arg.compression_level)?; 
        }
        Subcli::filter { read1, read2, nbase, length, complexity, average_qual, phred, chunk, thread, failed, out1, out2 } => {
            filter_fastq(&read1, &read2, nbase, length, complexity, average_qual, phred, chunk, thread, &failed, &out1, &out2, arg.compression_level)?;
        }
        Subcli::concat { read1, read2, out1, out2 } => {
            concat_fqstq_lane(&read1, &read2, &out1, &out2, arg.compression_level)?;
        }
        Subcli::remove { input, out, name , save, rm} => {
            remove_read(input.as_ref(), out.as_ref(), &name, &save, rm, arg.compression_level)?;
        }
        Subcli::rename { input, keep, prefix, label, before, output } => {
            rename_fastq(input.as_ref(), keep, prefix, label.as_ref(), before, output.as_ref(), arg.compression_level)?;
        }
        Subcli::reverse { input, rev, out } => {
            reverse_comp_seq(input.as_ref(), out.as_ref(), rev, arg.compression_level)?;
        }
        Subcli::split { input, pre,   out, } => {
            split_interleaved(input.as_ref(), &out, &pre, arg.compression_level)?;
        }
        Subcli::merge { read1,   read2,  out, } => {
            interleaved(&read1, &read2, out.as_ref(), arg.compression_level)?;
        }
        Subcli::mask { input, phred, low, chars, out } => {
            mask_fastq(input.as_ref(), phred, low, chars,  out.as_ref(), arg.compression_level)?;
        }
        Subcli::split2 { input, num, gzip, bzip2, xz, name } => {
            split_chunk(input.as_ref(), num, gzip, bzip2, xz, &name, arg.compression_level)?;
        }
        Subcli::gcplot { input, output, show, prefix, width, height, ylim, types } => {
            gc_content(input.as_ref(), output.as_ref(), show, prefix, width, height, ylim, &types, arg.compression_level)?;
        }
        Subcli::length { input, out } => {
            fq_length(input.as_ref(), out.as_ref(), arg.compression_level)?;
        }
        Subcli::view { input, out } => {
            view_fq(input.as_ref(), out.as_ref(), arg.compression_level)?;
        }
    }

    Ok(())
}