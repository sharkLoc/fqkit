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
fqkit: a simple program for fastq file manipulation

Usage: fqkit <COMMAND>

Commands:
  subfq
          subsample sequences from big fastq file
  search
          search reads/motifs from fastq file
  stats
          summary for fastq format file
  plot
          line plot for fastq quaility stats
  fq2fa
          translate fastq to fasta
  barcode
          split barcode for PE reads
  remove
          remove reads by read name
  split
          split interleaved fastq file
  merge
          merge PE reads as interleaved fastq file
  gcplot
          get GC content result and plot
  help
          Print this message or the help of the given subcommand(s)

Options:
  -h, --help
          Print help information
  -V, --version
          Print version information

```
