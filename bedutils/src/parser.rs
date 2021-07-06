use sequencetoolkit_common::{SequenceToolkitError, SequenceToolkitErrorKind};
use std::io;
use std::str;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BedRegion {
    pub chromosome: Vec<u8>,
    pub start: u64,
    pub end: u64,
    pub columns: Vec<Vec<u8>>,
}

impl BedRegion {
    pub fn new() -> Self {
        BedRegion::default()
    }
}

pub struct BedWriter<R: io::Write> {
    writer: R,
}

impl<R: io::Write> BedWriter<R> {
    pub fn new(writer: R) -> Self {
        BedWriter { writer }
    }

    pub fn write_record(&mut self, region: &BedRegion) -> Result<(), SequenceToolkitError> {
        self.writer.write_all(&region.chromosome)?;
        write!(self.writer, "\t{}\t{}", region.start, region.end)?;
        for one in &region.columns {
            self.writer.write_all(b"\t")?;
            self.writer.write_all(one)?;
        }
        self.writer.write_all(b"\n")?;
        Ok(())
    }
}

pub struct BedReader<R: io::BufRead> {
    reader: R,
    buffer: Vec<u8>,
    line: u64,
}

impl<R: io::BufRead> BedReader<R> {
    pub fn new(reader: R) -> Self {
        BedReader {
            reader,
            buffer: Vec::new(),
            line: 0,
        }
    }

    pub fn next(&mut self, region: &mut BedRegion) -> Result<bool, SequenceToolkitError> {
        loop {
            self.line += 1;
            self.buffer.clear();
            self.reader.read_until(b'\n', &mut self.buffer)?;
            if self.buffer.is_empty() {
                return Ok(false);
            }
            if !self.buffer.starts_with(b"#") {
                break;
            }
        }

        if self.buffer.ends_with(b"\r\n") {
            self.buffer.pop();
            self.buffer.pop();
        } else if self.buffer.ends_with(b"\n") {
            self.buffer.pop();
        }

        let mut max = 0;
        for (i, v) in self.buffer.split(|x| *x == b'\t').enumerate() {
            match i {
                0 => {
                    region.chromosome.clear();
                    region.chromosome.extend_from_slice(&v);
                    region.start = 0;
                    region.end = 0;
                }
                1 => {
                    region.start = str::from_utf8(&v)?.parse()?;
                }
                2 => {
                    region.end = str::from_utf8(&v)?.parse()?;
                }
                _ => {
                    if region.columns.len() <= i - 3 {
                        region.columns.push(v.to_vec());
                    } else {
                        region.columns[i - 3].clear();
                        region.columns[i - 3].extend_from_slice(&v);
                    }
                }
            }
            max = i;
        }

        if max < 2 {
            return Err(SequenceToolkitErrorKind::BedParseError(self.line).into());
        }

        if max < region.columns.len() + 2 {
            region.columns.drain((max - 2)..);
        }

        Ok(true)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_bed_parser() -> Result<(), SequenceToolkitError> {
        let test_bed = b"X\t1\t2\tA\r\nY\t30\t40\nM\t500\t600\ta\tb\tc";
        let mut bed_reader = BedReader::new(&test_bed[..]);
        let mut region = BedRegion::default();
        bed_reader.next(&mut region)?;
        assert_eq!(
            region,
            BedRegion {
                chromosome: b"X".to_vec(),
                start: 1,
                end: 2,
                columns: vec![b"A".to_vec()]
            }
        );
        bed_reader.next(&mut region)?;
        assert_eq!(
            region,
            BedRegion {
                chromosome: b"Y".to_vec(),
                start: 30,
                end: 40,
                columns: vec![]
            }
        );
        bed_reader.next(&mut region)?;
        assert_eq!(
            region,
            BedRegion {
                chromosome: b"M".to_vec(),
                start: 500,
                end: 600,
                columns: vec![b"a".to_vec(), b"b".to_vec(), b"c".to_vec(),]
            }
        );
        Ok(())
    }
}
