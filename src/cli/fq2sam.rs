use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::info;
use paraseq::fastq;

#[allow(clippy::too_many_arguments)]
pub fn fastq2sam(
    r1: &String,
    r2: Option<&String>,
    sm: &str,
    rg: Option<String>,
    lb: Option<String>,
    pl: Option<String>,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    info!("sample name set: {}", sm);

    let mut sam = file_writer(out, compression_level, stdout_type)?;
    let rg = if let Some(x) = rg {
        x
    } else {
        String::from("A")
    };
    sam.write_fmt(format_args!(
        "@HD\tVN:1.6\tSO:queryname\n@RG\tID:{}\tSM:{}",
        rg, sm
    ))?;
    if let Some(lb) = lb {
        sam.write_fmt(format_args!("\tLB:{}", lb))?;
        info!("library name set: {}", lb);
    }
    if let Some(pl) = pl {
        sam.write_fmt(format_args!("\tPL:{}", pl))?;
        info!("platform set: {}", pl);
    }
    sam.write_all("\n".as_bytes())?;

    if let Some(r2) = r2 {
        let mut fq1 = file_reader(Some(r1)).map(fastq::Reader::new)?;
        let mut fq2 = file_reader(Some(r2)).map(fastq::Reader::new)?;

        let mut rset1 = fastq::RecordSet::default();
        let mut rset2 = fastq::RecordSet::default();

        while rset1.fill(&mut fq1)? && rset2.fill(&mut fq2)? {
            for (rec1, rec2) in rset1
                .iter()
                .map_while(Result::ok)
                .zip(rset2.iter().map_while(Result::ok))
            {
                sam.write_all(rec1.id())?;
                sam.write_all(b"\t77\t*\t0\t0\t*\t*\t0\t0\t")?;
                sam.write_all(rec1.sep())?;
                sam.write_all(b"\t")?;
                sam.write_all(rec1.qual())?;
                sam.write_all(b"\tRG:Z:")?;
                sam.write_all(rg.as_bytes())?;
                sam.write_all(b"\n")?;

                sam.write_all(rec2.id())?;
                sam.write_all(b"\t141\t*\t0\t0\t*\t*\t0\t0\t")?;
                sam.write_all(rec2.sep())?;
                sam.write_all(b"\t")?;
                sam.write_all(rec2.qual())?;
                sam.write_all(b"\tRG:Z:")?;
                sam.write_all(rg.as_bytes())?;
                sam.write_all(b"\n")?;
            }
        }
    } else {
        let mut fq = file_reader(Some(r1)).map(fastq::Reader::new)?;

        let mut rset = fastq::RecordSet::default();
        while rset.fill(&mut fq)? {
            for rec in rset.iter().map_while(Result::ok) {
                sam.write_all(rec.id())?;
                sam.write_all(b"\t4\t*\t0\t0\t*\t*\t0\t0\t")?;
                sam.write_all(rec.sep())?;
                sam.write_all(b"\t")?;
                sam.write_all(rec.qual())?;
                sam.write_all(b"\tRG:Z:")?;
                sam.write_all(rg.as_bytes())?;
                sam.write_all(b"\n")?;
            }
        }
    }
    sam.flush()?;
    Ok(())
}
