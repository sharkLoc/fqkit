use clap::{Parser,value_parser};

#[derive(Parser, Debug)]
#[command(
    name = "FqKit",
    author = "sharkLoc",
    version = "0.3.11",
    about = "A simple and cross-platform program for fastq file manipulation",
    long_about = None,
    next_line_help = false,
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
#[command(help_template =
    "{name} -- {about}\n\nVersion: {version}\
    \n\nAuthors: {author} <mmtinfo@163.com>\
    \nSource code: https://github.com/sharkLoc/fqkit.git\
    \n\n{before-help}
{usage-heading} {usage}\n\n{all-args}\n\nUse \"fqkit help [command]\" for more information about a command"
)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Subcli,
    /// set gzip/bzip2/xz compression level 1 (compress faster) - 9 (compress better) for gzip/bzip2/xz output file,
    /// just work with option -o/--out
    #[arg(long = "compress-level", default_value_t = 6, global = true,
        value_parser = value_parser!(u32).range(1..=9), value_name = "int",
        help_heading = Some("Global Arguments")
    )]
    pub compression_level: u32,
    /// if file name specified, write log message to this file, or write to stderr
    #[arg(long = "log", global = true, help_heading = Some("Global Arguments"), value_name = "str")]
    pub logfile: Option<String>,
    /// control verbosity of logging, possible values: {error, warn, info, debug, trace}
    #[arg(short = 'v', long = "verbosity", global = true, value_name = "str",
        default_value_t = String::from("debug"),
        help_heading = Some("Global Arguments")
    )]
    pub verbose: String,
    /// be quiet and do not show any extra information
    #[arg(short = 'q', long = "quiet", global= true, help_heading = Some("Global FLAGS"))]
    pub quiet: bool,
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
        #[arg(short = 'n', long = "num", default_value_t = 10)]
        num: usize,
        /// output fastq file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// get last N records from fastq file
    tail {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// print first N fastq records
        #[arg(short = 'n', long = "num", default_value_t = 10)]
        num: usize,
        /// output fastq file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// concat fastq files from different lanes
    concat {
        /// input read1 list file, one fastq file per line
        #[arg(short = 'i', long = "input1")]
        read1: String,
        /// input read2 list file, one fastq file per line
        #[arg(short = 'I', long = "input2")]
        read2: String,
        /// read1 output file name,  files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out1")]
        out1: String,
        /// read2 output file name,  files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'O', long = "out2")]
        out2: String,
    },
    /// subsample sequences from big fastq file.
    #[command(visible_alias = "sample")]
    subfq {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// set rand seed.
        #[arg(short = 's', long = "seed", default_value_t = 69)]
        seed: u64,
        /// subseq number
        #[arg(short = 'n', long = "num")]
        num: usize,
        /// read files twice to reduce much memory but cost more time
        #[arg(short = 'r', long = "rdc", help_heading = Some("FLAGS"))]
        rdc: bool,
        /// fastq output file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// select pair-end reads by read id
    select {
        /// input read1 fastq file
        #[arg(short = '1', long = "read1")]
        read1: String,
        /// input read2 fastq file
        #[arg(short = '2', long = "read2")]
        read2: String,
        /// output selected  forward(read1) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short='f', long = "out1")]
        out1: String,
        /// output selected resverse(read2) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short='r', long = "out2")]
        out2: String,
    },
    /// trim fastq reads by position
    trim {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// trim int bp from left  
        #[arg(short, long, default_value_t=0)]
        left: usize,
        /// trim int bp from right
        #[arg(short, long, default_value_t=0)]
        right: usize,
        /// fastq output file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// a simple filter for pair end fastq sqeuence
    filter {
        /// input read1 fastq file
        #[arg(short = '1', long = "read1")]
        read1: String,
        /// input read2 fastq file
        #[arg(short = '2', long = "read2")]
        read2: String,
        /// if one read number of N base is more then N base limit, then this read pair is discarded.
        #[arg(short = 'n', long = "n-limit", default_value_t=5)]
        nbase: usize,
        /// reads shorter than length_required will be discarded
        #[arg(short = 'l', long = "length", default_value_t=30)]
        length: usize,
        /// the complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]),
        ///a 51-bp sequence, with 3 bases that is different from its next base
        ///seq = 'AAAATTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGCCCC',  and complexity = 3/(51-1) = 6%
        ///the threshold for low complexity filter (0~100). 30 is recommended, which means 30% complexity is required.
        #[arg(short = 'y', long = "complexity", default_value_t = 0,
            value_parser = value_parser!(u32).range(0..=100),
            verbatim_doc_comment
        )]
        complexity: u32,
        /// if one read's average quality score < average qual, then this read pair is discarded,
        ///eg. Q20 error 0.01, Q30 error 0.001, averaging the probability of error is 0.0055 => Q value 22.59637
        #[arg(short = 'q', long = "average_qual", default_value_t = 20, verbatim_doc_comment)]
        average_qual: u8,
        ///phred score 33 or 64
        #[arg(short = 'p', long = "phred", default_value_t = 33)]
        phred: u8,
        /// the number of reads in the chunk on each thread
        #[arg(short, long, default_value_t = 5000)]
        chunk: usize,
        /// number of additional worker threads to use
        #[arg(short='@', long="thread", default_value_t = 4)]
        thread: usize,
        /// specify the file to store reads(interleaved) that cannot pass the filters, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short='u', long = "failed")]
        failed: String,
        /// output pass filtered  forward(read1) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short='f', long = "out1")]
        out1: String,
        /// output pass filtered resverse(read2) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short='r', long = "out2")]
        out2: String,
    },
    /// print fastq records in a range
    range {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// skip first int read records
        #[arg(short = 's', long = "skip", default_value_t = 0)]
        skip: usize,
        /// take int read records
        #[arg(short = 't', long = "take")]
        take: usize,
        /// fastq output file name or write to stdout, files ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,

    },
    /// search reads/motifs from fastq file
    search {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// specify pattern/motif, regular expression supported, e.g., -p "ATC{2,}" or -p "ATCCG"
        ///for multiple motifs, -p "TTAGGG|CCCTAA"
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
        /// output contain pattern/motif reads result fastq file or write to stdout,
        ///file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", verbatim_doc_comment)]
        out: Option<String>,
    },
    /// grep fastq sequence by read id or full name
    grep {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// read name list file, one name per line and without read name prefix "@"
        #[arg(short = 'i', long = "id-list")]
        ids: String,
        /// if specified, match read by full name instead of just id
        #[arg(short = 'f', long = "--full-name", help_heading = Some("FLAGS"))]
        full: bool,
        /// output matched reads result in fastq file or write to stdout,
        ///file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out", verbatim_doc_comment)]
        out: Option<String>,
    },
    /// summary for fastq format file
    #[command(visible_alias = "stat", subcommand_help_heading = Some("Statistics"))]
    stats {
        /// input fastq file, or read from stdin
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
    /// shuffle fastq sequences 
    #[command(before_help = "note: all records will be readed into memory")]
    shuffle {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// set rand seed.
        #[arg(short = 's', long = "seed", default_value_t = 69)]
        seed: u64,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// report the number sequences and bases
    size {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// the number of reads in the chunk on each thread
        #[arg(short, long, default_value_t = 5000)]
        chunk: usize,
        /// number of additional worker threads to use
        #[arg(short='@', long="thread", default_value_t = 3)]
        thread: usize,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,

    },
    /// extract subsequences in sliding windows
    slide {
        /// input fastq file, or read from stdin
        input: Option<String>,
        ///set sliding window step size
        #[arg(short = 'w', long = "window", default_value_t = 10)]
        window: usize,
        ///set sliding window step size
        #[arg(short = 's', long = "step", default_value_t = 5)]
        step: usize,
        /// suffix added to the sequence ID
        #[arg(short = 'S', long = "suffidx", default_value_t = String::from("_slide"))]
        suffix: String,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
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
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// if specified, remove sequence id description
        #[arg(short='r', long="remove", help_heading = Some("FLAGS"))]
        remove: bool,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// converts a fastq file to an unaligned SAM file
    fq2sam {
        /// input fastq file
        #[arg(short = '1', long = "read1")]
        r1: String,
        /// input fastq file for the second read of paired end data
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
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
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
        #[arg(short = 'o', long = "out")]
        out: Option<String>,

    },
    /// flatten fastq sequences
    #[command(visible_alias = "flat")]
    flatten {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// filed number, id:1, sequence:2, symbol:4, quality:8
        ///eg. output id, sequence and quality value: 1 + 2 + 8 == 11 ,
        #[arg(short = 'f', long = "field", default_value_t = 3, verbatim_doc_comment)]
        flag: u8,
        /// output seprater, can be ",",  ";", 
        #[arg(short = 's', long = "sep", default_value_t='\t')]
        sep: char,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// join paired end reads that are overlapping into a single longer read
    #[command(before_help=r"Note:
    1. if overlapping regions are detected, low quality bases are corrected by high quality paired bases. 
    2. if a base is corrected, the quality value is also corrected.
    3. only paired end reads such as the following will be detected and merged, overlap mode:
       r1: GCAAGCGTTAATCGGAATTTATGGGCGTAAAGCGCACGCAGGA
                            |\|||||||||||||||||||||||\ 
                            TCTATGGGCGTAAAGCGCACGCAGGCATGCTGGGCGTAAAGCGCACGCAGGC  r2: reverse complement
    ")]
    join {
        /// input read1 fastq file
        #[arg(short = '1', long = "read1")]
        read1: String,
        /// input read2 fastq file
        #[arg(short = '2', long = "read2")]
        read2: String,
        /// minimum overlap length in PE reads
        #[arg(short = 'l', long = "length", default_value_t=15)]
        length: usize,
        /// maximum mismatch count in overlap region
        #[arg(short = 'm', long = "miss", default_value_t=10)]
        miss: usize,
        /// the number of reads in the chunk on each thread
        #[arg(short, long, default_value_t = 5000)]
        chunk: usize,
        /// number of additional worker threads to use
        #[arg(short='@', long="thread", default_value_t = 6)]
        thread: usize,
        /// output unmerged forward(read1) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short='f', long = "out1")]
        out1: String,
        /// output unmerged resverse(read2) fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short='r', long = "out2")]
        out2: String,
        /// output merged fastq file name,  file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "merged")]
        merged: String,
    },
    /// perform demultiplex for pair-end fastq reads 
    barcode {
        /// input read1 fastq file
        #[arg(short = '1', long = "read1")]
        read1: String,
        /// input read2 fastq file <barcode in this file>
        #[arg(short = '2', long = "read2")]
        read2: String,
        /// barcode list file, format eg:
        ///ATGCAGTG    sample1
        ///TGCAGTAC    sample2
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
        #[arg(short = 'o', long = "outdir", default_value_t = String::from("."))]
        outdir: String,
    },
    /// check the validity of a fastq record
    #[command(before_help = "note: this function will return an Err if one of the following conditions is met:\n
      1. the record identifier is empty.
      2. there is a non-ASCII character found in either the sequence or quality strings.
      3. the sequence and quality strings do not have the same length.\n")]
    check {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// if set, just save correct reads
        #[arg(short = 's', long = "save", help_heading = Some("FLAGS"))]
        save: bool, 
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// remove reads by read name.
    remove {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
        /// read name list file, one name per line and without read name prefix "@"
        #[arg(short = 'n', long = "name")]
        name: String,
        /// save removed reads in read name list
        #[arg(short = 's', long = "save",default_value_t=String::from("rm.fq.gz"))]
        save: String,
    },
    /// rename sequence id in fastq file
    rename {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// if specified, keep sequence id description
        #[arg(short = 'k', long = "keep", help_heading = Some("FLAGS"))]
        keep: bool, 
        /// set new id prefix for sequence
        #[arg(short = 'p', long = "prefix")]
        prefix: Option<String>,
        /// output fastq file name, or write to stdout, file name ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
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
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// split interleaved fastq file
    split {
        /// input fastq file, or read from stdin
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
        /// input read1 fastq file.
        #[arg(short = '1', long = "read1")]
        read1: String,
        /// input read2 fastq file.
        #[arg(short = '2', long = "read2")]
        read2: String,
        /// output interleaved fastq file name, eg. result.fq.bz2
        #[arg(short = 'o', long = "out", default_value_t = String::from("interleaved.fq.gz"))]
        out: String,
    },
    /// convert any low quality base to 'N' or other chars 
    mask {
        /// input fastq file, or read from stdin
        input: Option<String>,
        ///phred score 33 or 64
        #[arg(short = 'p', long = "phred", default_value_t = 33)]
        phred: u8,
        /// low quality
        #[arg(short = 'l', long = "low-quality",default_value_t = 5)] 
        low: u8,
        /// mask low quality ( <= low quality) base with this char
        #[arg(short = 'c', long = "char", default_value_t = 'N')]
        chars: char,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// split fastq file by records number
    split2 {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// set record number for each mini fastq file
        #[arg(short = 'n', long = "num", default_value_t = 200000)]
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
        /// output prefix name
        #[arg(short = 'p', long = "prefix", default_value_t = String::from("sub"))]
        name: String,
    },
    /// get GC content result and plot
    gcplot {
        /// input fastq file, or read from stdin
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
    /// get reads length count
    length {
        /// input fastq file, or read from stdin
        input: Option<String>,
        /// output file name or write to stdout, file ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    },
    /// view fastq file page by page
    view {
        /// input fastq file
        input: Option<String>,
        /// output reads page by page, file name ending in .gz/.bz2/.xz will be compressed automatically
        #[arg(short = 'o', long = "out")]
        out: Option<String>,
    }
}

