use clap::{
    ArgAction, Parser,
    builder::{
        Styles,
        styling::{AnsiColor, Effects},
    },
    value_parser,
};

// Configures Clap v3-style help menu colors
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default());

#[derive(Parser, Debug)]
#[command(styles = STYLES)]
#[command(
    name = "FqKit",
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = "A simple and cross-platform program for fastq file manipulation",
    long_about = None,
    next_line_help = false,
    disable_help_flag = true,
    disable_version_flag = true,
    propagate_version = true,
    before_help = r"Fqkit supports reading and writing gzip (.gz) format.
Bzip2 (.bz2) format is supported since v0.3.8.
Xz (.xz) format is supported since v0.3.9.
Under the same compression level, xz has the highest compression ratio but consumes more time. 

Compression level:
  format   range   default   crate
  gzip     1-9     6         https://crates.io/crates/flate2
  bzip2    1-9     6         https://crates.io/crates/bzip2
  xz       1-9     6         https://crates.io/crates/xz2"
)]
#[command(help_template = "{name} -- {about}\n\nVersion: {version}\
    \n\nAuthors: {author} <mmtinfo@163.com>\
    \nSource code: https://github.com/sharkLoc/fqkit.git\
    \n\n{before-help}
{usage-heading} {usage}\n\n{all-args}\n\nUse \"fqkit help [command]\" for more information about a command")]
pub struct Args {
    #[clap(subcommand)]
    pub command: Subcli,

    /// threads number
    #[arg(short = '@', long = "threads", default_value_t = 4, global = true, value_name = "INT", help_heading = Some("Global Arguments"))]
    pub threads: usize,

    /// set gzip/bzip2/xz compression level 1 (compress faster) - 9 (compress better) for gzip/bzip2/xz output file, just work with option -o/--out
    #[arg(long = "compress-level", default_value_t = 6, global = true,
        value_parser = value_parser!(u32).range(1..=9), value_name = "INT", help_heading = Some("Global Arguments")
    )]
    pub compression_level: u32,

    /// output type for stdout: 'g' gzip; 'b' bzip2; 'x' xz; 'u' uncompressed txt format
    #[arg(long = "output-type", global = true, help_heading = Some("Global Arguments"), value_name = "u|g|b|x", default_value_t = 'u')]
    pub stdout_type: char,

    /// if file name specified, write log message to this file, or write to stderr
    #[arg(long = "log", global = true, help_heading = Some("Global Arguments"), value_name = "FILE")]
    pub logfile: Option<String>,

    /// control verbosity of logging, [-v: Error, -vv: Warn, -vvv: Info, -vvvv: Debug, -vvvvv: Trace, defalut: Debug]
    #[arg(short = 'v', long = "verbosity", action = ArgAction::Count, global = true,
        default_value_t = 4, help_heading = Some("Global Arguments")
    )]
    pub verbose: u8,

    /// be quiet and do not show any extra information
    #[arg(short = 'q', long = "quiet", global= true, help_heading = Some("Global FLAGS"))]
    pub quiet: bool,

    /// prints help information
    #[arg(short = 'h', long, action = ArgAction::Help, global= true, help_heading = Some("Global FLAGS"))]
    pub help: Option<String>,

    /// prints version information
    #[arg(short = 'V', long, action = ArgAction::Version, global= true, help_heading = Some("Global FLAGS"))]
    pub version: Option<String>,
}

#[derive(Parser, Debug)]
#[allow(non_camel_case_types)]
pub enum Subcli {
    /// get first N records from fastq file
    #[command(visible_alias = "head")]
    topn {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// print first N fastq records
        #[arg(short = 'n', long = "num", default_value_t = 10, value_name = "INT")]
        num: usize,
        /// output fastq file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// get last N records from fastq file
    #[command(
        before_help = "note: read files twice to reduce much memory but cost more time, can't use in Stream"
    )]
    tail {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// print last N fastq records
        #[arg(short = 'n', long = "num", default_value_t = 10, value_name = "INT")]
        num: usize,
        /// output fastq file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// concat fastq files from different lanes
    concat {
        /// input read1 list file, one fastq file per line
        #[arg(short = 'i', long = "input1", value_name = "FILE")]
        read1: String,
        /// input read2 list file, one fastq file per line
        #[arg(short = 'I', long = "input2", value_name = "FILE")]
        read2: String,
        /// read1 output file name,  files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out1", value_name = "FILE")]
        out1: String,
        /// read2 output file name,  files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'O', long = "out2", value_name = "FILE")]
        out2: String,
    },
    /// subsample sequences from big fastq file.
    #[command(visible_alias = "sample")]
    subfq {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// set rand seed.
        #[arg(short = 's', long = "seed", default_value_t = 69, value_name = "INT")]
        seed: u64,
        /// subseq number
        #[arg(short = 'n', long = "num", value_name = "INT")]
        num: usize,
        /// read files twice to reduce much memory but cost more time
        #[arg(short = 'r', long = "rdc", help_heading = Some("FLAGS"))]
        rdc: bool,
        /// fastq output file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// select pair-end reads by read id
    select {
        /// input read1 fastq file
        #[arg(short = '1', long = "read1", value_name = "FILE")]
        read1: String,
        /// input read2 fastq file
        #[arg(short = '2', long = "read2", value_name = "FILE")]
        read2: String,
        /// output selected  forward(read1) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'f', long = "out1", value_name = "FILE")]
        out1: String,
        /// output selected resverse(read2) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'r', long = "out2", value_name = "FILE")]
        out2: String,
    },
    /// trim fastq reads by position
    trim {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// trim int bp from left  
        #[arg(short, long, default_value_t = 0, value_name = "INT")]
        left: usize,
        /// trim int bp from right
        #[arg(short, long, default_value_t = 0, value_name = "INT")]
        right: usize,
        /// after trimming, reads shorter than INT are discarded
        #[arg(short = 'd', long = "discard", default_value_t = 0, value_name = "INT")]
        len: usize,
        /// fastq output file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// cut the adapter sequence on the reads
    adapter {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// adapter sequence file in Fasta format
        #[arg(short = 'f', long = "fasta", value_name = "FILE")]
        fa: String,
        /// addpter on the left side of the read
        #[arg(short, long, help_heading = Some("FLAGS"))]
        left: bool,
        /// max missmatch allowed
        #[arg(short, long, default_value_t = 0, value_name = "INT")]
        miss: usize,
        /// fastq output file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// a simple filter for pair end fastq sqeuence
    filter {
        /// input read1 fastq file
        #[arg(short = '1', long = "read1", value_name = "FILE")]
        read1: String,
        /// input read2 fastq file
        #[arg(short = '2', long = "read2", value_name = "FILE")]
        read2: String,
        /// if one read number of N base is more then N base limit, then this read pair is discarded.
        #[arg(short = 'n', long = "n-limit", default_value_t = 5, value_name = "INT")]
        nbase: usize,
        /// reads shorter than length required will be discarded
        #[arg(short = 'l', long = "length", default_value_t = 30, value_name = "INT")]
        length: usize,
        /// the complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]), a 51-bp sequence, with 3 bases that is different from its next base
        ///seq = 'AAAATTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGCCCC' and complexity = 3/(51-1) = 6%, the threshold for low complexity filter (0~100). 30 is recommended, which means 30% complexity is required.
        #[arg(short = 'y', long = "complexity", default_value_t = 0, value_parser = value_parser!(u32).range(0..=100), value_name = "INT")]
        complexity: u32,
        /// if one read's average quality score < average qual, then this read pair is discarded,
        ///eg. Q20 error 0.01, Q30 error 0.001, averaging the probability of error is 0.0055 => Q value 22.59637
        #[arg(
            short = 'Q',
            long = "average_qual",
            default_value_t = 20,
            value_name = "INT"
        )]
        average_qual: u8,
        ///phred score 33 or 64
        #[arg(short = 'p', long = "phred", default_value_t = 33, value_name = "INT")]
        phred: u8,
        /// if set, specify the file to store reads(interleaved) that cannot pass the filters, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'u', long = "failed", value_name = "FILE")]
        failed: Option<String>,
        /// output pass filtered  forward(read1) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'f', long = "out1", value_name = "FILE")]
        out1: String,
        /// output pass filtered resverse(read2) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'r', long = "out2", value_name = "FILE")]
        out2: String,
    },
    /// join paired end reads that are overlapping into a single longer read
    #[command(before_help = r"Note:
    1. if overlapping regions are detected, low quality bases are corrected by high quality paired bases. 
    2. if a base is corrected, the quality value is also corrected.
    3. only paired end reads such as the following will be detected and merged, overlap mode:
       r1: GCAAGCGTTAATCGGAATTTATGGGCGTAAAGCGCACGCAGGA
                            |\|||||||||||||||||||||||\ 
                            TCTATGGGCGTAAAGCGCACGCAGGCATGCTGGGCGTAAAGCGCACGCAGGC  r2: reverse complement                       
    ")]
    join {
        /// input read1 fastq file
        #[arg(short = '1', long = "read1", value_name = "FILE")]
        read1: String,
        /// input read2 fastq file
        #[arg(short = '2', long = "read2", value_name = "FILE")]
        read2: String,
        /// minimum overlap length in PE reads
        #[arg(short = 'l', long = "length", default_value_t = 30)]
        length: usize,
        /// maximum mismatch rate count in overlap region
        #[arg(short = 'm', long = "miss", default_value_t = 0.1)]
        miss: f64,
        /// output joinde long fastq file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "output")]
        output: Option<String>,
        /// output interleaved fastq file name for non-overlap pe reads
        #[arg(short = 'n', long = "non-overlap", value_name = "FILE")]
        non: Option<String>,
    },
    /// print fastq records in a range
    range {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// skip first int read records
        #[arg(short = 's', long = "skip", default_value_t = 0, value_name = "INT")]
        skip: usize,
        /// take int read records
        #[arg(short = 't', long = "take", value_name = "INT")]
        take: usize,
        /// fastq output file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// search reads/motifs from fastq file
    search {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// specify pattern/motif, regular expression supported, e.g., -p "ATC{2,}" or -p "ATCCG", for multiple motifs, -p "TTAGGG|CCCTAA"
        #[arg(short = 'p', long = "pattern", value_name = "STR")]
        pat: String,
        /// if specified, enable case insensitive matching for the entire pattern
        #[arg(short = 'i', long = "ignore-case", help_heading = Some("FLAGS"))]
        case: bool,
        /// invert the sense of matching, to select non-matching reads
        #[arg(short = 'u', long = "invert-match", help_heading = Some("FLAGS"))]
        invert: bool,
        /// output contain pattern/motif reads result fastq file or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// grep fastq sequence by read id or full name
    grep {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// read name list file, one name per line and without read name prefix "@"
        #[arg(short = 'i', long = "id-list", value_name = "FILE")]
        ids: String,
        /// if specified, match read by full name instead of just id
        #[arg(short = 'f', long = "full-name", help_heading = Some("FLAGS"))]
        full: bool,
        /// output matched reads result in fastq file or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// summary for fastq format file
    #[command(visible_alias = "stat", subcommand_help_heading = Some("Statistics"))]
    stats {
        /// input fastq file, or read from stdin
        input: Option<String>,
        ///phred score 33 or 64
        #[arg(short = 'p', long = "phred", default_value_t = 33, value_name = "INT")]
        phred: u8,
        /// specify a name for summary output file
        #[arg(short='s',long="sumy",default_value_t=String::from("summary.txt"), value_name = "FILE")]
        sum: String,
        /// if not specified, cycle result write to stdout
        #[arg(short = 'c', long = "cycle", value_name = "FILE")]
        cyc: Option<String>,
    },
    /// a simple kmer counter
    kmer {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// set kmer size
        #[arg(
            short = 'k',
            long = "kmer-size",
            default_value_t = 21,
            value_name = "INT"
        )]
        size: usize,
        /// add header info in output file
        #[arg(short = 'H', long, help_heading = Some("FLAGS"))]
        header: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// shuffle fastq sequences
    #[command(before_help = "note: all records will be readed into memory")]
    shuffle {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// set rand seed.
        #[arg(short = 's', long = "seed", default_value_t = 69, value_name = "INT")]
        seed: u64,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// report the number sequences and bases
    size {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// extract subsequences in sliding windows
    slide {
        /// input fastq file, or read from stdin
        input: Option<String>,
        ///set sliding window step size
        #[arg(short = 'w', long = "window", default_value_t = 10, value_name = "INT")]
        window: usize,
        ///set sliding window step size
        #[arg(short = 's', long = "step", default_value_t = 5, value_name = "INT")]
        step: usize,
        /// suffix added to the sequence ID
        #[arg(short = 'S', long = "suffidx", default_value_t = String::from("_slide"), value_name = "STR")]
        suffix: String,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// sort fastq file by name/seq/gc/length
    #[command(before_help = "note: all records will be readed into memory")]
    sort {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// sort reads by name
        #[arg(short = 'n', long = "sort-by-name" ,help_heading = Some("FLAGS"))]
        name: bool,
        /// sort reads by sequence
        #[arg(short = 's', long = "sort-by-seq" ,help_heading = Some("FLAGS"))]
        seq: bool,
        /// sort reads by gc content
        #[arg(short = 'g', long = "sort-by-gc", help_heading = Some("FLAGS"))]
        gc: bool,
        /// sort read by length
        #[arg(short = 'l', long = "sort-by-length", help_heading = Some("FLAGS"))]
        length: bool,
        /// output reversed result
        #[arg(short = 'r', long = "reverse", help_heading = Some("FLAGS"))]
        reverse: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// line plot for A T G C N percentage in read position
    plot {
        /// input cycle result data: fqkit stats cycle output
        #[arg(short = 'd', long = "data", value_name = "FILE")]
        data: String,
        /// if specified, show line plot in terminal
        #[arg( short = 's', long ="show-terminal", help_heading = Some("FLAGS"))]
        show: bool,
        /// output base figure prefix name
        #[arg(short='p', long="prefix", default_value_t=String::from("base_plot"), value_name = "STR")]
        prefix: String,
        /// set output figure width
        #[arg(short = 'W', long = "width", default_value_t = 960, value_name = "INT")]
        width: usize,
        /// set output figure height
        #[arg(
            short = 'H',
            long = "height",
            default_value_t = 540,
            value_name = "INT"
        )]
        height: usize,
        /// set max ylim (0~100)
        #[arg(
            short = 'y',
            long = "ylim",
            default_value_t = 50.0,
            value_name = "FLOAT"
        )]
        ylim: f32,
        /// figure type 'png' or 'svg'
        #[arg(short='t', long="types", default_value_t=String::from("png"), value_name = "STR")]
        types: String,
    },
    /// translate fastq to fasta
    fq2fa {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// if specified, remove sequence id description
        #[arg(short='r', long="remove", help_heading = Some("FLAGS"))]
        remove: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// converts a fastq file to an unaligned SAM file
    fq2sam {
        /// input fastq file
        #[arg(short = '1', long = "read1", value_name = "FILE")]
        r1: String,
        /// input fastq file for the second read of paired end data
        #[arg(short = '2', long = "read2", help_heading = Some("Optional Arguments"), value_name = "FILE")]
        r2: Option<String>,
        /// sample name to insert into the read group header
        #[arg(short = 's', long = "sample-name", value_name = "STR")]
        sm: String,
        /// read group name, default: A
        #[arg(short = 'r', long = "read-group-name", help_heading = Some("Optional Arguments") ,value_name = "STR")]
        rg: Option<String>,
        /// the library name to place into the LB attribute in the read group header
        #[arg(short = 'l', long = "library-name", help_heading = Some("Optional Arguments") ,value_name = "STR")]
        lb: Option<String>,
        /// the platform type (e.g. ILLUMINA, SOLID) to insert into the read group header
        #[arg(short = 'p', long = "platform", help_heading = Some("Optional Arguments") ,value_name = "STR")]
        pl: Option<String>,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// converts the fastq file quality scores
    fqscore {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// converts the quality scores from phred 64 to phred 33, quality - 31
        #[arg(long = "to33", help_heading = Some("FLAGS"))]
        to33: bool,
        /// converts the quality scores from phred 33 to phred 64, quality + 31
        #[arg(long = "to64", help_heading = Some("FLAGS"))]
        to64: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// flatten fastq sequences
    #[command(visible_alias = "flat")]
    flatten {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// filed number, id:1, sequence:2, symbol:4, quality:8; eg. output id, sequence and quality value: 1 + 2 + 8 == 11
        #[arg(short = 'f', long = "field", default_value_t = 3, value_name = "INT")]
        flag: u8,
        /// output seprater, can be ",",  ";",
        #[arg(short = 's', long = "sep", default_value_t = '\t', value_name = "CHAR")]
        sep: char,
        /// if specified, add N base count in output
        #[arg(short = 'n', long = "gap-n", help_heading = Some("FLAGS"))]
        gap: bool,
        /// if specified, add read length in output
        #[arg(short = 'l', long = "length", help_heading = Some("FLAGS"))]
        len: bool,
        /// if specified, add GC content(%) in output
        #[arg(short = 'g', long = "gc-content", help_heading = Some("FLAGS"))]
        gc: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// perform demultiplex for pair-end fastq reads
    #[command(visible_alias = "demux")]
    barcode {
        /// input read1 fastq file
        #[arg(short = '1', long = "read1", value_name = "FILE")]
        read1: String,
        /// input read2 fastq file <barcode in this file>
        #[arg(short = '2', long = "read2", value_name = "FILE")]
        read2: String,
        /// barcode list file, format eg:
        /// ATGCAGTG    sample1
        /// TGCAGTAC    sample2
        #[arg(
            short = 'b',
            long = "barcode",
            verbatim_doc_comment,
            value_name = "FILE"
        )]
        bar: String,
        /// barcode position mode, 1:left, 2:right
        #[arg(short = 'm', long = "mode", default_value_t = 2, value_name = "INT")]
        mode: usize,
        /// barcode reverse complement
        #[arg(short = 'r', long = "rev_comp", help_heading = Some("FLAGS"))]
        trans: bool,
        /// barcode mismatch base count
        #[arg(short = 'e', long = "error", default_value_t = 0, value_name = "INT")]
        mismatch: usize,
        /// if specified, output gzip compressed file
        #[arg(short = 'z', long = "gzip", help_heading = Some("FLAGS"))]
        gzip: bool,
        /// if specified, output bzip2 compressed file
        #[arg(short = 'Z', long = "bzip2", help_heading = Some("FLAGS"))]
        bzip2: bool,
        /// if specified, output xz compressed file
        #[arg(short = 'x', long = "xz", help_heading = Some("FLAGS"))]
        xz: bool,
        /// fastq file output dir.
        #[arg(short = 'o', long = "outdir", default_value_t = String::from(".") ,value_name = "DIR")]
        outdir: String,
    },
    /// remove reads by read name.
    #[command(visible_alias = "rm")]
    remove {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// read name list file, one name per line and without read name prefix "@"
        #[arg(short = 'n', long = "name", value_name = "FILE")]
        name: String,
        /// save removed reads in read name list
        #[arg(short = 's', long = "save",default_value_t=String::from("rm.fq.gz") ,value_name = "FILE")]
        save: String,
        /// if set, do not output removed reads
        #[arg(short = 'r', long = "remove", help_heading = Some("FLAGS"))]
        rm: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// rename sequence id in fastq file
    #[command(visible_alias = "rn")]
    rename {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// if specified, keep sequence id description
        #[arg(short = 'k', long = "keep", help_heading = Some("FLAGS"))]
        keep: bool,
        /// set new id prefix for sequence
        #[arg(short = 'p', long = "prefix", value_name = "STR")]
        prefix: Option<String>,
        /// if specified, add a label before/after read id
        #[arg(short = 'l', long = "label", value_name = "STR")]
        label: Option<String>,
        /// if set, label before read id
        #[arg(short = 'b', long = "before", help_heading = Some("FLAGS"))]
        before: bool,
        /// output fastq file name, or write to stdout, file name ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "STR")]
        output: Option<String>,
    },
    /// get a reverse-complement of fastq file.
    #[command(visible_alias = "rev")]
    reverse {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// if set, just output reverse sequences, the quality scores are also reversed
        #[arg(short = 'r', long = "reverse", help_heading = Some("FLAGS"))]
        rev: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// split interleaved fastq file
    split {
        /// input fastq file, or read from stdin
        //#[arg(short = 'i', long = "input" ,value_name = "STR")]
        input: Option<String>,
        /// if specified, output gzip compressed file
        #[arg(short = 'z', long = "gzip", help_heading = Some("FLAGS"))]
        gzip: bool,
        /// if specified, output bzip2 compressed file
        #[arg(short = 'Z', long = "bzip2", help_heading = Some("FLAGS"))]
        bzip2: bool,
        /// if specified, output xz compressed file
        #[arg(short = 'x', long = "xz", help_heading = Some("FLAGS"))]
        xz: bool,
        /// output fastq file prefix name
        #[arg(short = 'p', long = "prefix" , default_value_t = String::from("demo"), value_name = "STR")]
        pre: String,
        /// fastq file outdir
        #[arg(short = 'o', long = "outdir", default_value_t = String::from(".") ,value_name = "DIR")]
        out: String,
    },
    /// merge PE reads as interleaved fastq file
    merge {
        /// input read1 fastq file.
        #[arg(short = '1', long = "read1", value_name = "FILE")]
        read1: String,
        /// input read2 fastq file.
        #[arg(short = '2', long = "read2", value_name = "FILE")]
        read2: String,
        /// output interleaved fastq file name, eg. result.fq.bz2
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// convert any low quality base to 'N' or other chars
    mask {
        /// input fastq file, or read from stdin
        input: Option<String>,
        ///phred score 33 or 64
        #[arg(short = 'p', long = "phred", default_value_t = 33, value_name = "INT")]
        phred: u8,
        /// low quality
        #[arg(
            short = 'l',
            long = "low-quality",
            default_value_t = 5,
            value_name = "INT"
        )]
        low: u8,
        /// mask low quality ( <= low quality) base with this char
        #[arg(short = 'c', long = "char", default_value_t = 'N', value_name = "CHAR")]
        chars: char,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// split fastq file by records number
    split2 {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// set record number for each mini fastq file
        #[arg(
            short = 'n',
            long = "num",
            default_value_t = 200000,
            value_name = "INT"
        )]
        num: usize,
        /// if specified, output gzip compressed file
        #[arg(short = 'z', long = "gzip", help_heading = Some("FLAGS"))]
        gzip: bool,
        /// if specified, output bzip2 compressed file
        #[arg(short = 'Z', long = "bzip2", help_heading = Some("FLAGS"))]
        bzip2: bool,
        /// if specified, output xz compressed file
        #[arg(short = 'x', long = "xz", help_heading = Some("FLAGS"))]
        xz: bool,
        /// output file prefix name
        #[arg(short = 'p', long = "prefix", default_value_t = String::from("sub") ,value_name = "STR")]
        name: String,
        /// fastq file output dir.
        #[arg(short = 'o', long = "outdir", default_value_t = String::from(".") ,value_name = "DIR")]
        outdir: String,
    },
    /// get GC content result and plot
    #[command(visible_alias = "gc")]
    gcplot {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// output GC contnet result file name
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        output: Option<String>,
        /// if specified, show histogram graphs in terminal
        #[arg( short = 's', long ="show-terminal", help_heading = Some("FLAGS"))]
        show: bool,
        /// output base figure prefix name
        #[arg(short='p', long="prefix", default_value_t=String::from("gc_plot") ,value_name = "STR")]
        prefix: String,
        /// set output figure width
        #[arg(short = 'W', long = "width", default_value_t = 960, value_name = "INT")]
        width: usize,
        /// set output figure height
        #[arg(
            short = 'H',
            long = "height",
            default_value_t = 540,
            value_name = "INT"
        )]
        height: usize,
        /// set max ylim (0~100)
        #[arg(short = 'y', long = "ylim", default_value_t = 15, value_name = "INT")]
        ylim: usize,
        /// figure type 'png' or 'svg'
        #[arg(short='t', long="types", default_value_t=String::from("png") ,value_name = "STR")]
        types: String,
    },
    /// get reads length count
    #[command(visible_alias = "len")]
    length {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// output reversed result
        #[arg(short = 'r', long = "reverse", help_heading = Some("FLAGS"))]
        reverse: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
    /// view fastq file page by page
    view {
        /// input fastq file
        input: Option<String>,
        /// output reads page by page, file name ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", value_name = "FILE")]
        out: Option<String>,
    },
}
