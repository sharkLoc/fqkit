use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use colored::*;
use log::error;
use std::io;
use term_size::dimensions;

pub fn view_fq(
    file: Option<&String>,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    if file.is_none() {
        error!("do not read file from stdin.");
        std::process::exit(1);
    }

    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut fq_writer = file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;

    let mut iter_fq = fq_reader.records().map_while(Result::ok).peekable();
    let (mut page, mut start, mut end) = (0usize, 1usize, 0usize);

    loop {
        let term_height = dimensions().unwrap().1;
        let page_size = term_height - 2;
        let this = page_size / 4;
        for _ in 0..this {
            if let Some(rec) = iter_fq.next() {
                fq_writer.write_record(&rec)?;
                fq_writer.flush()?;
            } else {
                eprintln!("{}", "End of file!".red());
                return Ok(());
            }
        }
        page += 1;
        end += this;
        eprint!("{}", format!("Page: {page}, reads from {start} to {end} , press enter to show next page, q and enter to quit: ").green());
        start += this;

        let mut input = String::new();
        io::stdin().read_line(&mut input)?;

        match input.trim() {
            "q" => {
                break;
            }
            _ => {
                continue;
            }
        }
    }

    Ok(())
}
