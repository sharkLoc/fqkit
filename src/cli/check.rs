use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::{error, info};

pub fn check_fastq(
    file: Option<&String>,
    save: bool,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let (mut total, mut ok_read, mut err_read) = (0, 0, 0);
    let fp_reader = file_reader(file).map(fastq::Reader::new)?;

    if save {
        let mut out_writer =
            file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;
        for rec in fp_reader.records().map_while(Result::ok) {
            total += 1;
            match rec.check() {
                Ok(_) => {
                    ok_read += 1;
                    out_writer.write_record(&rec)?;
                }
                Err(s) => {
                    err_read += 1;
                    error!("record order: {}, id: {}, {}", total, rec.id(), s);
                }
            };
        }
        out_writer.flush()?;
    } else {
        for rec in fp_reader.records().map_while(Result::ok) {
            total += 1;
            match rec.check() {
                Ok(_) => {
                    ok_read += 1;
                }
                Err(s) => {
                    err_read += 1;
                    error!("record order: {}, id: {}, {}", total, rec.id(), s);
                }
            };
        }
    }

    info!(
        "total reads num: {}, ok reads number: {}, error reads number: {}",
        total, ok_read, err_read
    );

    Ok(())
}
