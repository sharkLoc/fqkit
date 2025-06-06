use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::{fastq, fastq::Record};
use log::{info, trace};

pub fn slide_fastq(
    file: Option<&String>,
    step: usize,
    wind: usize,
    out: Option<&String>,
    suffix: &str,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let fq_reader = file_reader(file).map(fastq::Reader::new)?;
    info!("window size : {}", wind);
    info!("step size: {}", step);
    let mut fq_writer = file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;
    let mut window = wind;
    for rec in fq_reader.records().map_while(Result::ok) {
        let seq = rec.seq();
        let qual = rec.qual();
        let len = seq.len();
        let mut st = 0;
        loop {
            if window < len {
                let this_desc = if let Some(desc) = rec.desc() {
                    format!("{}{}:{}-{}", desc, suffix, st + 1, window)
                } else {
                    format!("{}:{}-{}", suffix, st + 1, window)
                };
                fq_writer.write_record(&Record::with_attrs(
                    rec.id(),
                    Some(this_desc.as_str()),
                    &seq[st..window],
                    &qual[st..window],
                ))?;
                st += step;
                window += step;
            } else {
                if st < len {
                    let this_desc = if let Some(desc) = rec.desc() {
                        format!("{}{}:{}-{}", desc, suffix, st + 1, len)
                    } else {
                        format!("{}:{}-{}", suffix, st + 1, len)
                    };
                    fq_writer.write_record(&Record::with_attrs(
                        rec.id(),
                        Some(this_desc.as_str()),
                        &seq[st..len],
                        &qual[st..len],
                    ))?;
                } else {
                    trace!(
                        "slice read start position: {} is bigger than read length: {}",
                        st + 1,
                        len
                    );
                }
                // single read slice done, init window size with input value
                window = wind;
                break;
            }
        }
    }
    fq_writer.flush()?;

    Ok(())
}
