use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::*;
use paraseq::{fastq, fastx::Record};

#[allow(clippy::too_many_arguments)]
pub fn rename_fastq(
    input: Option<&String>,
    keep: bool,
    prefix: Option<String>,
    label: Option<&String>,
    before: bool,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let mut fq_reader = fastq::Reader::new(file_reader(input)?);
    let mut rset = fastq::RecordSet::default();

    let mut writer = file_writer(output, compression_level, stdout_type)?;
    let mut n: usize = 0;

    while rset.fill(&mut fq_reader)? {
        for rec in rset.iter().map_while(Result::ok) {
            n += 1;

            let mut newid: Vec<Vec<u8>> = vec![];
            match &prefix {
                Some(pre) => {
                    let mut id_split = rec.id_str().split_whitespace();

                    if before {
                        let lab = label.unwrap_or(&String::new()).as_bytes().to_vec();
                        let pre = pre.as_bytes().to_vec();
                        let n = n.to_string().as_bytes().to_vec();
                        newid.push(lab);
                        newid.push(pre);
                        newid.push(n);

                        if keep {
                            id_split.next();

                            // add desc in new id
                            if id_split.clone().count() > 0 {
                                newid.push(" ".as_bytes().to_vec());
                                let desc = id_split
                                    .collect::<Vec<&str>>()
                                    .join(" ")
                                    .as_bytes()
                                    .to_vec();
                                newid.push(desc);
                            }
                        }
                    } else {
                        let lab = label.unwrap_or(&String::new()).as_bytes().to_vec();
                        let pre = pre.as_bytes().to_vec();
                        let n = n.to_string().as_bytes().to_vec();
                        newid.push(pre);
                        newid.push(lab);
                        newid.push(n);

                        if keep {
                            id_split.next();

                            // add desc in new id
                            if id_split.clone().count() > 0 {
                                newid.push(" ".as_bytes().to_vec());
                                let desc = id_split
                                    .collect::<Vec<&str>>()
                                    .join(" ")
                                    .as_bytes()
                                    .to_vec();
                                newid.push(desc);
                            }
                        }
                    }
                }
                None => {
                    let lab = label.unwrap_or(&String::new()).as_bytes().to_vec();
                    if before {
                        newid.push(lab);
                        if keep {
                            newid.push(rec.id().to_vec());
                        } else {
                            let id = rec
                                .id_str()
                                .split_whitespace()
                                .next()
                                .unwrap()
                                .as_bytes()
                                .to_vec();
                            newid.push(id);
                        }
                    } else {
                        let mut id_split = rec.id_str().split_whitespace();
                        if keep {
                            newid.push(id_split.next().unwrap().as_bytes().to_vec());
                            newid.push(lab);

                            // add desc in new id
                            if id_split.clone().count() > 0 {
                                newid.push(" ".as_bytes().to_vec());
                                let desc = id_split
                                    .collect::<Vec<&str>>()
                                    .join(" ")
                                    .as_bytes()
                                    .to_vec();
                                newid.push(desc);
                            }
                        } else {
                            newid.push(id_split.next().unwrap().as_bytes().to_vec());
                            newid.push(lab);
                        }
                    }
                }
            }

            let id = newid.concat();
            write_record(&mut writer, id.as_slice(), rec.seq(), rec.qual())?;
        }
    }
    writer.flush()?;

    info!("total rename sequence number: {}", n);
    Ok(())
}
