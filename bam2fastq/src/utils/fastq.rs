use super::{DNBSeqFASTQInfo, IlluminaFASTQInfo, DNBSEQ_REGEX_IN_BAM, ILLUMINA_REGEX_IN_BAM};
use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::io::{self, prelude::*};
use std::str;

pub type TilePosition = u32;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Deserialize, Serialize)]
pub struct GenericFastqEntry {
    read_name: Vec<u8>,
    sequence: Vec<u8>,
    quality: Vec<u8>,
}

impl GenericFastqEntry {
    pub fn read<R: BufRead>(mut reader: R) -> anyhow::Result<Option<GenericFastqEntry>> {
        let mut buffer = [0u8; 2];
        if reader.read(&mut buffer[0..1])? == 0 {
            return Ok(None);
        }
        if buffer[0] != b'@' {
            return Err(anyhow::anyhow!(
                "Invalid FASTQ name: {}",
                str::from_utf8(&buffer).unwrap()
            ));
        }

        let mut read_name = Vec::new();
        if reader.read_until(b'\n', &mut read_name)? == 0 {
            return Err(anyhow::anyhow!("unexpected EOF"));
        }
        read_name.remove(read_name.len() - 1);

        // load sequence
        let mut sequence = Vec::new();
        if reader.read_until(b'\n', &mut sequence)? == 0 {
            return Err(anyhow::anyhow!(
                "Uncompleted FASTQ 1: read name: {}",
                str::from_utf8(&read_name).unwrap()
            ));
        }
        sequence.remove(sequence.len() - 1);

        // load read name 2
        reader.read_exact(&mut buffer[0..2])?;
        if &buffer[0..2] != b"+\n" {
            return Err(anyhow::anyhow!(
                "Invalid FASTQ read name2: {} for read name: {}",
                str::from_utf8(&buffer).unwrap(),
                str::from_utf8(&read_name).unwrap()
            ));
        }

        // load quality
        let mut quality = Vec::new();
        if reader.read_until(b'\n', &mut quality)? == 0 {
            return Err(anyhow::anyhow!(
                "Uncompleted FASTQ 3: read name: {}",
                str::from_utf8(&read_name).unwrap()
            ));
        }
        quality.remove(quality.len() - 1);

        Ok(Some(GenericFastqEntry {
            read_name,
            sequence,
            quality,
        }))
    }

    pub fn write<W: Write>(&self, mut writer: W) -> io::Result<()> {
        writer.write_all(b"@")?;
        writer.write_all(&self.read_name)?;
        writer.write_all(b"\n")?;
        writer.write_all(&self.sequence)?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(&self.quality)?;
        writer.write_all(b"\n")?;
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Deserialize, Serialize)]
pub struct IlluminaFastqEntry {
    prefix: String,
    tile: u16,
    y_pos: TilePosition,
    x_pos: TilePosition,
    sequence: Vec<u8>,
    quality: Vec<u8>,
}

impl TryFrom<GenericFastqEntry> for IlluminaFastqEntry {
    type Error = anyhow::Error;
    fn try_from(entry: GenericFastqEntry) -> anyhow::Result<IlluminaFastqEntry> {
        let read_name_str = str::from_utf8(&entry.read_name).unwrap();
        if let Some(cap) = ILLUMINA_REGEX_IN_BAM.captures(read_name_str) {
            let prefix = cap.name("prefix").unwrap().as_str().to_string();
            let tile: u16 = cap.name("tile").unwrap().as_str().parse().unwrap();
            let x_pos: TilePosition = cap.name("x_pos").unwrap().as_str().parse().unwrap();
            let y_pos: TilePosition = cap.name("y_pos").unwrap().as_str().parse().unwrap();
            Ok(IlluminaFastqEntry {
                prefix,
                tile,
                x_pos,
                y_pos,
                sequence: entry.sequence,
                quality: entry.quality,
            })
        } else {
            Err(anyhow::anyhow!(
                "read name is not illumina format: {}",
                read_name_str
            ))
        }
    }
}

impl IlluminaFastqEntry {
    pub fn write<W: Write>(
        &self,
        mut writer: W,
        read_number: Option<&str>,
        fastq_info: Option<&IlluminaFASTQInfo>,
    ) -> io::Result<()> {
        write!(
            writer,
            "@{}:{}:{}:{}",
            self.prefix, self.tile, self.x_pos, self.y_pos
        )?;

        if let Some(fastq_info) = fastq_info {
            if self.prefix != fastq_info.prefix {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "FASTQ prefix is not match",
                ));
            }

            if let Some(special_index) = fastq_info
                .special_index
                .get(&(self.tile, self.x_pos, self.y_pos))
            {
                log::trace!(
                    "special index found for {:?}",
                    (self.tile, self.x_pos, self.y_pos)
                );
                write!(writer, " {}:{}", read_number.unwrap_or("?"), special_index)?;
            } else {
                log::trace!(
                    "No special index for {:?}",
                    (self.tile, self.x_pos, self.y_pos)
                );
                write!(
                    writer,
                    " {}:{}",
                    read_number.unwrap_or("?"),
                    fastq_info.common_index
                )?;
            }
        }
        writer.write_all(b"\n")?;
        writer.write_all(&self.sequence)?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(&self.quality)?;
        writer.write_all(b"\n")?;
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Deserialize, Serialize)]
pub struct DNBSeqFastqEntry {
    group_index: usize,
    index: String,
    sequence: Vec<u8>,
    quality: Vec<u8>,
}

impl DNBSeqFastqEntry {
    pub fn from_fastq_entry(
        entry: GenericFastqEntry,
        fastq_info: &DNBSeqFASTQInfo,
    ) -> anyhow::Result<Self> {
        let read_name_str = str::from_utf8(&entry.read_name).unwrap();
        if let Some(cap) = DNBSEQ_REGEX_IN_BAM.captures(read_name_str) {
            let prefix = cap.name("prefix").unwrap().as_str();
            if prefix != fastq_info.prefix {
                return Err(anyhow::anyhow!("Unknown prefix: {}", read_name_str));
            }
            let group_name = cap.name("group").unwrap().as_str();
            if let Some(group_index) = fastq_info.group_order_rev.get(group_name) {
                Ok(DNBSeqFastqEntry {
                    group_index: *group_index,
                    index: cap.name("index").unwrap().as_str().to_string(),
                    sequence: entry.sequence,
                    quality: entry.quality,
                })
            } else {
                Err(anyhow::anyhow!(
                    "group name is not found in FASTQ info: {}",
                    read_name_str
                ))
            }
        } else {
            Err(anyhow::anyhow!(
                "read name is not DNBSeq format: {}",
                read_name_str
            ))
        }
    }

    pub fn write<W: Write>(
        &self,
        mut writer: W,
        read_number: Option<&str>,
        fastq_info: &DNBSeqFASTQInfo,
    ) -> io::Result<()> {
        write!(
            writer,
            "@{}{}{}",
            fastq_info.prefix, fastq_info.group_order[self.group_index], self.index
        )?;

        if let Some(read_number) = read_number {
            write!(writer, "/{}", read_number)?;
        }

        writer.write_all(b"\n")?;
        writer.write_all(&self.sequence)?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(&self.quality)?;
        writer.write_all(b"\n")?;
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;
    #[test]
    fn test_lookup() {
        let mut map = HashMap::<(u16, u16, u16), String>::new();
        map.insert((1, 2, 3), "hoge".to_string());
        assert_eq!(map.get(&(1, 2, 3)).unwrap(), "hoge");
    }
}
