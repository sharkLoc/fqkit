# fqkit
🦀 a simple program for fastq file manipulation


## install
##### setp1：install cargo first 
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

##### step2:
```bash
cargo install fqkit
# or

git clone https://github.com/sharkLoc/fqkit.git
cd fqkit
cargo b --release
# mv target/release/fqkit to anywhere you want 
```

## usage

```bash
fqkit: A simple program for fastq file manipulation

Version: 0.2.20
Authors: sharkLoc <mmtinfo@163.com>

Usage: fqkit [OPTIONS] <COMMAND>

Commands:
  topn     get first N records from fastq file
  subfq    subsample sequences from big fastq file [aliases: sample]
  trim     trim fastq file
  range    print fastq records in a range
  search   search reads/motifs from fastq file
  stats    summary for fastq format file [aliases: stat]
  shuffle  shuffle fastq sequences
  size     report the number sequences and bases
  sort     sort fastq file by name/seq/gc/length
  plot     line plot for A T G C N percentage in read position
  fq2fa    translate fastq to fasta
  fq2sam   converts a fastq file to an unaligned SAM file
  flatten  flatten fastq sequences
  barcode  split barcode for PE reads
  check    check the validity of a fastq record
  remove   remove reads by read name
  reverse  get a reverse-complement of fastq file [aliases: rev]
  split    split interleaved fastq file
  merge    merge PE reads as interleaved fastq file
  mask     convert any low quality base to 'N' or other chars
  split2   split fastq file by records number
  gcplot   get GC content result and plot
  view     view fastq file page by page
  help     Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

Global Arguments:
  -v, --verbosity <VERBOSE>  control verbosity of logging, Possible values: {error, warn, info, debug, trace} [default: debug]
      --log <LOGFILE>        if specified, write log message to file, or write to stderr

Global FLAGS:
  -q, --quiet  be quiet and do not show extra information
```

## change log
<details>
<summary>timeline</summary>

2023.11.03:
 - update to version 0.2.12
 - add subcommand trim
 - update cmd help information

2023.11.08:
 - update to version 0.2.13
 - add subcommand reverse

2023.11.10:
 - update to version 0.2.14
 - add subcommand view
 - rebuilt some command interface

2023.11.29:
 - update to version 0.2.15
 - recode func in stats subcommand

2023.12.04
 - update to version 0.2.16
 - add subcommand size

2023.12.05
 - update to version 0.2.17
 - update code for subcommand size and search
 - add subcommand fq2sam
 - add log information options verbosity

2023.12.09
 - update to version 0.2.18
 - add subcommand sort, range, check and mask
 - update to version 0.2.19
 - add subcommand shuffle
</details>

#### ** any bugs please report issues **💖