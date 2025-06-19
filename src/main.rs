use anyhow::Error;
use clap::Parser;
use log::{debug, error, info};
use std::time::Instant;

mod logger;
use logger::*;
mod command;
mod errors;
mod utils;
use command::*;
mod cli;
use cli::{
    barcode::*, concat::*, cutadapter::*, filter::*, flatten::*, fq2fa::*, fq2sam::*, fqscore::*,
    gcplot::*, grep::*, join::*, kmer::*, length::*, mask::*, merge::*, plot::*, range::*,
    remove::*, rename::*, reverse::*, search::*, select::*, shuffle::*, size::*, slide::*, sort::*,
    split::*, split2::*, stats::*, subfq::*, tail::*, top::*, trimfq::*, view::*,
};

fn main() {
    match run_main() {
        Ok(_) => {}
        Err(e) => {
            error!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

fn run_main() -> Result<(), Error> {
    let arg = Args::parse();
    logger(arg.verbose, arg.logfile, arg.quiet)?;
    let start = Instant::now();
    info!("version: {}", env!("CARGO_PKG_VERSION"));

    let cpus = num_cpus::get();
    debug!("cpu numbers: {}", cpus);
    rayon::ThreadPoolBuilder::new()
        .num_threads(arg.threads)
        .build_global()?;

    match arg.command {
        Subcli::topn { input, num, out } => {
            top_n_records(
                input.as_ref(),
                num,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::tail { input, num, out } => {
            tail_n_records(
                input.as_ref(),
                num,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::subfq {
            input,
            seed,
            num,
            rdc,
            out,
        } => {
            subset_fastq(
                rdc,
                input.as_ref(),
                num,
                seed,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::select {
            read1,
            read2,
            out1,
            out2,
        } => {
            select_pe_fastq(
                &read1,
                &read2,
                &out1,
                &out2,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::trim {
            input,
            left,
            right,
            len,
            out,
        } => {
            trim_fq(
                input.as_ref(),
                left,
                right,
                len,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::adapter {
            input,
            fa,
            left,
            miss,
            out,
        } => {
            cut_adapter(
                input.as_ref(),
                &fa,
                left,
                miss,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::range {
            input,
            skip,
            take,
            out,
        } => {
            range_fastq(
                input.as_ref(),
                skip,
                take,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::search {
            input,
            pat,
            case,
            invert,
            out,
        } => {
            search_fq(
                input.as_ref(),
                &pat,
                case,
                invert,
                out.as_ref(),
                arg.threads,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::grep {
            input,
            ids,
            full,
            out,
        } => {
            grep_fastq(
                input.as_ref(),
                &ids,
                full,
                arg.threads,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::fq2fa { input, remove, out } => {
            fq2fa(
                input.as_ref(),
                remove,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::fq2sam {
            r1,
            r2,
            sm,
            rg,
            lb,
            pl,
            out,
        } => {
            fastq2sam(
                &r1,
                r2.as_ref(),
                &sm,
                rg,
                lb,
                pl,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::fqscore {
            input,
            to33,
            to64,
            out,
        } => {
            phred_score(
                input.as_ref(),
                out.as_ref(),
                to33,
                to64,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::flatten {
            input,
            flag,
            sep,
            gap,
            len,
            gc,
            out,
        } => {
            flatten_fq(
                input.as_ref(),
                out.as_ref(),
                flag,
                sep,
                gap,
                len,
                gc,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::plot {
            data,
            show,
            prefix,
            width,
            height,
            ylim,
            types,
        } => {
            let df = cycle_data(Some(&data))?;
            let _x = plot_line(df, show, prefix, width, height, ylim, &types);
        }
        Subcli::stats {
            input,
            phred,
            sum,
            cyc,
        } => {
            stat_fq(
                input.as_ref(),
                &sum,
                cyc.as_ref(),
                phred,
                arg.threads,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::shuffle { input, seed, out } => {
            shuffle_fastq(
                input.as_ref(),
                seed,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::size { input, out } => {
            size_fastq(
                input.as_ref(),
                arg.threads,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::slide {
            input,
            window,
            step,
            suffix,
            out,
        } => {
            slide_fastq(
                input.as_ref(),
                step,
                window,
                out.as_ref(),
                &suffix,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::sort {
            input,
            name,
            seq,
            gc,
            length,
            reverse,
            out,
        } => {
            sort_fastq(
                input.as_ref(),
                name,
                seq,
                gc,
                length,
                reverse,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::barcode {
            read1,
            read2,
            bar,
            mode,
            trans,
            mismatch,
            gzip,
            bzip2,
            xz,
            outdir,
        } => {
            split_fq(
                &read1,
                &read2,
                &bar,
                trans,
                mode,
                mismatch,
                &outdir,
                gzip,
                bzip2,
                xz,
                arg.compression_level,
            )?;
        }
        Subcli::filter {
            read1,
            read2,
            nbase,
            length,
            complexity,
            average_qual,
            phred,
            failed,
            out1,
            out2,
        } => {
            filter_fastq(
                &read1,
                &read2,
                nbase,
                length,
                complexity,
                average_qual,
                phred,
                arg.threads,
                failed.as_ref(),
                &out1,
                &out2,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::concat {
            read1,
            read2,
            out1,
            out2,
        } => {
            concat_fqstq_lane(
                &read1,
                &read2,
                &out1,
                &out2,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::remove {
            input,
            name,
            save,
            rm,
            out,
        } => {
            remove_read(
                input.as_ref(),
                out.as_ref(),
                &name,
                &save,
                rm,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::rename {
            input,
            keep,
            prefix,
            label,
            before,
            output,
        } => {
            rename_fastq(
                input.as_ref(),
                keep,
                prefix,
                label.as_ref(),
                before,
                output.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::reverse { input, rev, out } => {
            reverse_comp_seq(
                input.as_ref(),
                out.as_ref(),
                rev,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::split {
            input,
            gzip,
            bzip2,
            xz,
            pre,
            out,
        } => {
            split_interleaved(
                input.as_ref(),
                &out,
                &pre,
                gzip,
                bzip2,
                xz,
                arg.compression_level,
            )?;
        }
        Subcli::merge { read1, read2, out } => {
            interleaved(
                &read1,
                &read2,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::mask {
            input,
            phred,
            low,
            chars,
            out,
        } => {
            mask_fastq(
                input.as_ref(),
                phred,
                low,
                chars,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::split2 {
            input,
            num,
            gzip,
            bzip2,
            xz,
            name,
            outdir,
        } => {
            split_chunk(
                input.as_ref(),
                num,
                gzip,
                bzip2,
                xz,
                &name,
                &outdir,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::gcplot {
            input,
            output,
            show,
            prefix,
            width,
            height,
            ylim,
            types,
        } => {
            gc_content(
                input.as_ref(),
                output.as_ref(),
                show,
                prefix,
                width,
                height,
                ylim,
                &types,
                arg.threads,
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::length {
            input,
            reverse,
            out,
        } => {
            fq_length(
                input.as_ref(),
                reverse,
                arg.threads,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::view { input, out } => {
            view_fq(
                input.as_ref(),
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::join {
            read1,
            read2,
            length,
            miss,
            output,
            non,
        } => {
            join_overlap(
                &read1,
                &read2,
                miss,
                length,
                output.as_ref(),
                non.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
        Subcli::kmer {
            input,
            size,
            header,
            out,
        } => {
            kmer_count(
                input.as_ref(),
                size,
                header,
                out.as_ref(),
                arg.compression_level,
                arg.stdout_type,
            )?;
        }
    }

    info!("time elapsed is: {:?}", start.elapsed());

    Ok(())
}
