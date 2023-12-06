use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    author = "sharkLoc",
    version = "0.2.17",
    about = "A simple program for fastq file manipulation",
    long_about = None,
    next_line_help = false,
)]
#[command(help_template = 
    "{name}: {about}\n\nVersion: {version}\
    \nAuthors: {author} <mmtinfo@163.com>\
    \n\n{usage-heading} {usage}\n\n{all-args}\n"
)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Subcli,
    /// control verbosity of logging, Possible values: {error, warn, info, debug, trace}
    #[arg(short = 'v', long = "verbosity", global = true, default_value_t = String::from("debug"), help_heading = Some("Global Arguments"))]
    pub verbose: String,
    /// be quiet and do not show extra information
    #[arg(short = 'q', long = "quiet", global= true, help_heading = Some("Global FLAGS"))]
    pub quiet: bool,
}

#[derive(Parser, Debug)]
#[allow(non_camel_case_types)]
pub enum Subcli {
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
        /// specify pattern/motif, regular expression supported, e.g., -p "ATC{2,}" or -p "ATCCG"
        /// for multiple motifs, -p "TTAGGG|CCCTAA"
        #[arg(short = 'p', long = "pattern",verbatim_doc_comment)]
        pat: String,
        /// if specified,  enable case insensitive matching for the entire pattern
        #[arg(short = 'i', long ="ignore-case", help_heading = Some("FLAGS"))]
        case: bool,
        /// the number of reads in the chunk on each thread
        #[arg(short, long, default_value_t = 5000)]
        chunk: usize,
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
        /// the number of reads in the chunk on each thread
        #[arg(short, long, default_value_t = 5000)]
        chunk: usize,
        /// number of additional worker threads to use
        #[arg(short='@', long="thread", default_value_t = 3)]
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
    /// converts a fastq file to an unaligned SAM file
    fq2sam {
        /// input fastq[.gz] file
        #[arg(short = '1', long = "read1")]
        r1: String,
        /// input fastq[.gz] file for the second read of paired end data
        #[arg(short = '2', long = "read2", help_heading = Some("Optional Arguments"))]
        r2: Option<String>,
        /// sample name to insert into the read group header
        #[arg(short = 's', long = "sample-name")]
        sm: String,
        /// read group name, default: A
        #[arg(short = 'r', long = "read-group-name", help_heading = Some("Optional Arguments"))]
        rg: Option<String>,
        /// the library name to place into the LB attribute in the read group header
        #[arg(short = 'l', long = "library-name", help_heading = Some("Optional Arguments"))]
        lb: Option<String>,
        /// the platform type (e.g. ILLUMINA, SOLID) to insert into the read group header
        #[arg(short = 'p', long = "platform", help_heading = Some("Optional Arguments"))]
        pl: Option<String>,
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

