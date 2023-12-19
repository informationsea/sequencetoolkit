pub mod fastq;
pub mod largereorder;
pub mod ora;
pub mod progress;

use fastq::TilePosition;
use once_cell::sync::Lazy;
use regex::Regex;
use rust_htslib::bam::record::Aux;
use std::collections::HashMap;
use std::fmt;
use std::io::{self, prelude::*};

pub static ILLUMINA_REGEX: Lazy<Regex> = Lazy::new(|| {
    regex::Regex::new(r"^@(?P<prefix>[\w\-]+:\w+:\w+:\d):(?P<tile>\d+):(?P<x_pos>\d+):(?P<y_pos>\d+) (?P<read>\d):(?P<read_index>[NY]:\d+:.*)\n$").unwrap()
});
pub static DNBSEQ_REGEX: Lazy<Regex> = Lazy::new(|| {
    regex::Regex::new(
        r"^@(?P<prefix>(E|V)\d+B?L\dC)(?P<group>\d{3}R\d{3})(?P<index>\d+)/(?P<read>[12])\n$",
    )
    .unwrap()
});

pub static ILLUMINA_REGEX_IN_BAM: Lazy<Regex> = Lazy::new(|| {
    regex::Regex::new(
        r"^(?P<prefix>[\w\-]+:\w+:\w+:\d):(?P<tile>\d+):(?P<x_pos>\d+):(?P<y_pos>\d+)( .*)?$",
    )
    .unwrap()
});
pub static DNBSEQ_REGEX_IN_BAM: Lazy<Regex> = Lazy::new(|| {
    regex::Regex::new(
        r"^(?P<prefix>(E|V)\d+B?L\dC)(?P<group>\d{3}R\d{3})(?P<index>\d+)(/(?P<read>[12]))?$",
    )
    .unwrap()
});

pub fn aux_to_string<'a>(aux: &Aux<'a>) -> Option<&'a [u8]> {
    match aux {
        Aux::HexByteArray(x) => Some(x.as_bytes()),
        Aux::String(x) => Some(x.as_bytes()),
        _ => None,
    }
}

pub fn aux_to_str<'a>(aux: &Aux<'a>) -> Option<&'a str> {
    match aux {
        Aux::HexByteArray(x) => Some(*x),
        Aux::String(x) => Some(*x),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum NameType {
    Illumina,
    DNBSeq,
    Empty,
}

impl NameType {
    pub fn suggest(read1_name: &str, read2_name: &str) -> anyhow::Result<NameType> {
        if let Some(cap1) = ILLUMINA_REGEX.captures(read1_name) {
            if let Some(cap2) = ILLUMINA_REGEX.captures(read2_name) {
                if cap1.name("prefix").unwrap().as_str() == cap2.name("prefix").unwrap().as_str()
                    && cap1.name("read_index").unwrap().as_str()
                        == cap2.name("read_index").unwrap().as_str()
                    && cap1.name("read").unwrap().as_str() == "1"
                    && cap2.name("read").unwrap().as_str() == "2"
                {
                    return Ok(NameType::Illumina);
                }
            }
            return Err(anyhow::anyhow!(
                "Unmatched FASTQ read name (illumina): {} / {}",
                read1_name.trim(),
                read2_name.trim()
            ));
        }

        if let Some(cap1) = DNBSEQ_REGEX.captures(read1_name) {
            if let Some(cap2) = DNBSEQ_REGEX.captures(read2_name) {
                if cap1.name("prefix").unwrap().as_str() == cap2.name("prefix").unwrap().as_str()
                    && cap1.name("group").unwrap().as_str() == cap2.name("group").unwrap().as_str()
                    && cap1.name("index").unwrap().as_str() == cap2.name("index").unwrap().as_str()
                    && cap1.name("read").unwrap().as_str() == "1"
                    && cap2.name("read").unwrap().as_str() == "2"
                {
                    return Ok(NameType::DNBSeq);
                }
            }
            return Err(anyhow::anyhow!(
                "Unmatched FASTQ read name (DNBSeq): {} / {}",
                read1_name.trim(),
                read2_name.trim()
            ));
        }

        Err(anyhow::anyhow!(
            "Unknown FASTQ read name: {} / {}",
            read1_name.trim(),
            read2_name.trim()
        ))
    }
}

impl fmt::Display for NameType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                NameType::Illumina => "Illumina",
                NameType::DNBSeq => "DNBSeq",
                NameType::Empty => "Empty",
            }
        )
    }
}

pub struct TSVReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> TSVReader<R> {
    pub fn new(reader: R) -> Self {
        TSVReader { reader }
    }

    pub fn next(&mut self) -> io::Result<Vec<String>> {
        let mut line = String::new();
        if self.reader.read_line(&mut line)? == 0 {
            Ok(vec![])
        } else {
            Ok(line
                .trim_end_matches('\n')
                .split('\t')
                .map(|x| x.to_string())
                .collect())
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct IlluminaFASTQInfo {
    pub fastq1: String,
    pub fastq2: String,
    pub fastq_type: String,
    pub prefix: String,
    pub common_index: String,
    pub special_index: HashMap<(u16, TilePosition, TilePosition), String>,
    pub fastq1_md5: String,
    pub fastq2_md5: String,
}

impl IlluminaFASTQInfo {
    pub fn load(reader: impl BufRead) -> anyhow::Result<Self> {
        let mut tsv_reader = TSVReader::new(reader);
        let mut fastq_info = Self::default();

        loop {
            let row = tsv_reader.next()?;
            if row.is_empty() {
                break;
            }
            match row[0].as_str() {
                "FASTQ1" => fastq_info.fastq1 = row[1].to_string(),
                "FASTQ2" => fastq_info.fastq2 = row[1].to_string(),
                "FASTQ_TYPE" => fastq_info.fastq_type = row[1].to_string(),
                "READ_PREFIX" => fastq_info.prefix = row[1].to_string(),
                "COMMON_INDEX" => fastq_info.common_index = row[1].to_string(),
                "FASTQ1_MD5" => fastq_info.fastq1_md5 = row[1].to_string(),
                "FASTQ2_MD5" => fastq_info.fastq2_md5 = row[1].to_string(),
                "SPECIAL_INDEX" => {
                    fastq_info.special_index.insert(
                        (
                            row[1].parse().unwrap(),
                            row[2].parse().unwrap(),
                            row[3].parse().unwrap(),
                        ),
                        row[4].to_string(),
                    );
                }
                "VERIFY_RESULT" => {
                    if row[1] != "OK" {
                        return Err(anyhow::anyhow!("Bad verify result: {:?}", row));
                    }
                }
                _ => return Err(anyhow::anyhow!("Unknown row: {:?}", row)),
            }
        }

        if fastq_info.fastq1.is_empty() {
            return Err(anyhow::anyhow!("FASTQ1 path is not found in info"));
        }
        if fastq_info.fastq2.is_empty() {
            return Err(anyhow::anyhow!("FASTQ2 path is not found in info"));
        }
        if fastq_info.prefix.is_empty() {
            return Err(anyhow::anyhow!("Read prefix is not found in info"));
        }
        if fastq_info.common_index.is_empty() {
            return Err(anyhow::anyhow!("Common index is not found in info"));
        }
        if fastq_info.fastq1_md5.is_empty() {
            return Err(anyhow::anyhow!("FASTQ1 MD5 is not found in info"));
        }
        if fastq_info.fastq2_md5.is_empty() {
            return Err(anyhow::anyhow!("FASTQ2 MD5 is not found in info"));
        }

        //log::debug!("fastq info: {:?}", fastq_info);
        Ok(fastq_info)
    }
}

#[derive(Debug, Clone, Default)]
pub struct DNBSeqFASTQInfo {
    pub fastq1: String,
    pub fastq2: String,
    pub fastq_type: String,
    pub prefix: String,
    pub group_order: Vec<String>,
    pub group_order_rev: HashMap<String, usize>,
    pub fastq1_md5: String,
    pub fastq2_md5: String,
}

impl DNBSeqFASTQInfo {
    pub fn load(reader: impl BufRead) -> anyhow::Result<Self> {
        let mut tsv_reader = TSVReader::new(reader);
        let mut fastq_info = Self::default();

        let mut group_order_index = 0;

        loop {
            let row = tsv_reader.next()?;
            if row.is_empty() {
                break;
            }
            match row[0].as_str() {
                "FASTQ1" => fastq_info.fastq1 = row[1].to_string(),
                "FASTQ2" => fastq_info.fastq2 = row[1].to_string(),
                "FASTQ_TYPE" => fastq_info.fastq_type = row[1].to_string(),
                "PREFIX" => fastq_info.prefix = row[1].to_string(),
                "GROUP_ORDER" => {
                    fastq_info.group_order.push(row[1].to_string());
                    fastq_info
                        .group_order_rev
                        .insert(row[1].to_string(), group_order_index);
                    group_order_index += 1;
                }
                "FASTQ1_MD5" => fastq_info.fastq1_md5 = row[1].to_string(),
                "FASTQ2_MD5" => fastq_info.fastq2_md5 = row[1].to_string(),
                "VERIFY_RESULT" => {
                    if row[1] != "OK" {
                        return Err(anyhow::anyhow!("Bad verify result: {:?}", row));
                    }
                }
                _ => return Err(anyhow::anyhow!("Unknown row: {:?}", row)),
            }
        }

        if fastq_info.fastq1.is_empty() {
            return Err(anyhow::anyhow!("FASTQ1 path is not found in info"));
        }
        if fastq_info.fastq2.is_empty() {
            return Err(anyhow::anyhow!("FASTQ2 path is not found in info"));
        }
        if fastq_info.prefix.is_empty() {
            return Err(anyhow::anyhow!("Read prefix is not found in info"));
        }
        if fastq_info.group_order.is_empty() {
            return Err(anyhow::anyhow!("Group order is not found in info"));
        }
        if fastq_info.fastq1_md5.is_empty() {
            return Err(anyhow::anyhow!("FASTQ1 MD5 is not found in info"));
        }
        if fastq_info.fastq2_md5.is_empty() {
            return Err(anyhow::anyhow!("FASTQ2 MD5 is not found in info"));
        }

        Ok(fastq_info)
    }
}

pub struct DigestReader<R: Read, D: digest::Digest> {
    reader: R,
    digest: D,
}

impl<R: Read, D: digest::Digest + digest::FixedOutputReset> DigestReader<R, D> {
    pub fn new(reader: R, digest: D) -> Self {
        DigestReader { reader, digest }
    }

    pub fn finalize_reset(&mut self) -> digest::Output<D> {
        self.digest.finalize_reset()
    }
}

impl<R: Read, D: digest::Digest> Read for DigestReader<R, D> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let result = self.reader.read(buf)?;
        self.digest.update(&buf[..result]);
        Ok(result)
    }
}

pub struct DigestWriter<W: Write, D: digest::Digest> {
    writer: W,
    digest: D,
}

impl<W: Write, D: digest::Digest> Write for DigestWriter<W, D> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let wrote_size = self.writer.write(buf)?;
        self.digest.update(&buf[..wrote_size]);
        Ok(wrote_size)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

impl<W: Write, D: digest::Digest + digest::FixedOutputReset> DigestWriter<W, D> {
    pub fn new(writer: W, digest: D) -> Self {
        DigestWriter { writer, digest }
    }

    pub fn finalize_reset(&mut self) -> digest::Output<D> {
        self.digest.finalize_reset()
    }
}
