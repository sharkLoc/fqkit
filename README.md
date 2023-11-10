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

Version: 0.2.14
Authors: sharkLoc <mmtinfo@163.com>

Usage: fqkit [OPTIONS] <COMMAND>

Commands:
  topn     get first N records from fastq file
  subfq    subsample sequences from big fastq file
  trim     trim fastq file
  search   search reads/motifs from fastq file
  stats    summary for fastq format file [aliases: stat]
  plot     line plot for A T G C N percentage in read position
  fq2fa    translate fastq to fasta
  flatten  flatten fastq sequences
  barcode  split barcode for PE reads
  remove   remove reads by read name
  reverse  get a reverse-complement of fastq file
  split    split interleaved fastq file
  merge    merge PE reads as interleaved fastq file
  split2   split fastq file by records number
  gcplot   get GC content result and plot
  view     view fastq file page by page
  help     Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help information
  -V, --version  Print version information

Global FLAGS:
  -q, --quiet  be quiet and do not show extra information
```

## change log
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