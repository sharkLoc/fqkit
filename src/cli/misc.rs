use std::io::{Result, Write};

pub fn write_record<W>(writer: &mut W, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()>
where
    W: Write + Send,
{
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => *b, // unknow base like N
        })
        .collect()
}
