use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{info, trace};
use paraseq::fastq;

pub fn slide_fastq(
    file: Option<&String>,
    step: usize,
    wind: usize,
    out: Option<&String>,
    suffix: &str,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut fq_reader = file_reader(file).map(fastq::Reader::new)?;
    let mut rset = fastq::RecordSet::default();

    info!("window size : {}", wind);
    info!("step size: {}", step);

    let mut writer = file_writer(out, compression_level, stdout_type)?;
    let mut window = wind;
    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            let seq = rec.seq();
            let qual = rec.qual();
            let len = seq.len();
            let mut st = 0;

            loop {
                if window < len {
                    let mut newid = vec![];
                    newid.push(rec.id().to_vec());
                    newid.push(
                        format!("{}:{}-{}", suffix, st + 1, window)
                            .as_bytes()
                            .to_vec(),
                    );
                    write_record(
                        &mut writer,
                        newid.concat().as_slice(),
                        &seq[st..window],
                        &qual[st..window],
                    )?;
                    st += step;
                    window += step;
                } else {
                    if st < len {
                        let mut newid = vec![];
                        newid.push(rec.id().to_vec());
                        newid.push(format!("{}:{}-{}", suffix, st + 1, len).as_bytes().to_vec());
                        write_record(
                            &mut writer,
                            newid.concat().as_slice(),
                            &seq[st..len],
                            &qual[st..len],
                        )?;
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
    }

    writer.flush()?;

    Ok(())
}
