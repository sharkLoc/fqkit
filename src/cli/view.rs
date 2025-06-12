use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use colored::*;
use log::error;
use paraseq::fastq;
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

    let mut fq_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut rset = fastq::RecordSet::default();
    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;

    while rset.fill(&mut fq_reader)? {
        let mut iter_fq = rset.iter().map_while(Result::ok).peekable();

        let (mut page, mut start, mut end) = (0usize, 1usize, 0usize);
        loop {
            let term_height = dimensions().unwrap().1;
            let page_size = term_height - 2;
            let this = page_size / 4;
            for _ in 0..this {
                if let Some(rec) = iter_fq.next() {
                    write_record(&mut fq_writer, rec.id(), rec.seq(), rec.qual())?;
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
    }
    Ok(())
}
