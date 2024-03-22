![icon](https://github.com/sharkLoc/fqkit/blob/main/doc/fqkit_icon.PNG)
<!-- ![icon](doc/fqkit_icon.PNG) -->

# fqkit
![Static Badge](https://img.shields.io/badge/Author-sharkLoc-blue)
![Static Badge](https://img.shields.io/badge/Tool-fqkit-red)
![Crates.io (latest)](https://img.shields.io/crates/dv/fqkit?labelColor=rgb&color=hex&link=https%3A%2F%2Fcrates.io%2Fcrates%2Ffqkit)
![Crates.io](https://img.shields.io/crates/d/fqkit?label=Total%20download%20in%20crate.io)
![GitHub Gist last commit](https://img.shields.io/github/gist/last-commit/a4910923a230b8975218a188528463d7?logo=github)



🦀 a simple program for fastq file manipulation


## install
##### setp1： install cargo first 
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

##### step2:  on linux or windows
```bash
cargo install fqkit
# or

git clone https://github.com/sharkLoc/fqkit.git
cd fqkit
cargo b --release
# mv target/release/fqkit to anywhere you want 
```
##### install latest version

```bash
cargo install --git https://github.com/sharkLoc/fqkit.git
```

## usage

```bash
FqKit -- A simple and cross-platform program for fastq file manipulation

Version: 0.3.12

Authors: sharkLoc <mmtinfo@163.com>
Source code: https://github.com/sharkLoc/fqkit.git

Fqkit supports reading and writing gzip (.gz) format.
Bzip2 (.bz2) format is supported since v0.3.8.
Xz (.xz) format is supported since v0.3.9.
Under the same compression level, xz has the highest compression ratio but consumes more time. 

Compression level:
  format   range   default   crate
  gzip     1-9     6         https://crates.io/crates/flate2
  bzip2    1-9     6         https://crates.io/crates/bzip2
  xz       1-9     6         https://crates.io/crates/xz2


Usage: fqkit [OPTIONS] <COMMAND>

Commands:
  topn     get first N records from fastq file [aliases: head]
  tail     get last N records from fastq file
  concat   concat fastq files from different lanes
  subfq    subsample sequences from big fastq file [aliases: sample]
  select   select pair-end reads by read id
  trim     trim fastq reads by position
  filter   a simple filter for pair end fastq sqeuence
  range    print fastq records in a range
  search   search reads/motifs from fastq file
  grep     grep fastq sequence by read id or full name
  stats    summary for fastq format file [aliases: stat]
  shuffle  shuffle fastq sequences
  size     report the number sequences and bases
  slide    extract subsequences in sliding windows
  sort     sort fastq file by name/seq/gc/length
  plot     line plot for A T G C N percentage in read position
  fq2fa    translate fastq to fasta
  fq2sam   converts a fastq file to an unaligned SAM file
  fqscore  converts the fastq file quality scores
  flatten  flatten fastq sequences [aliases: flat]
  join     join paired end reads that are overlapping into a single longer read
  barcode  perform demultiplex for pair-end fastq reads
  check    check the validity of a fastq record
  remove   remove reads by read name
  rename   rename sequence id in fastq file
  reverse  get a reverse-complement of fastq file [aliases: rev]
  split    split interleaved fastq file
  merge    merge PE reads as interleaved fastq file
  mask     convert any low quality base to 'N' or other chars
  split2   split fastq file by records number
  gcplot   get GC content result and plot
  length   get reads length count
  view     view fastq file page by page
  help     Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

Global Arguments:
      --compress-level <int>  set gzip/bzip2/xz compression level 1 (compress faster) - 9 (compress better) for gzip/bzip2/xz output file, just work with option -o/--out [default: 6]
      --log <str>             if file name specified, write log message to this file, or write to stderr
  -v, --verbosity <str>       control verbosity of logging, possible values: {error, warn, info, debug, trace} [default: debug]

Global FLAGS:
  -q, --quiet  be quiet and do not show any extra information

Use "fqkit help [command]" for more information about a command
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

2023.12.10
 - update to version 0.3.0
 - add glob option --compress-level for gzip output file

2023.12.11
- update to version 0.3.1
- add subcommand grep and fqscore

2023.12.11
- update to version 0.3.2
- add subcommand slide

2023.12.19
- update to version 0.3.5
- add subcommand filter
</details>

#### ** any bugs please report issues **💖
