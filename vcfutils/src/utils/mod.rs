pub mod recalc_af;
pub mod tablewriter;

use crate::error::{VCFUtilsError, VCFUtilsErrorKind};
use csv::Reader as CSVReader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use vcf::{U8Vec, VCFReader};

type AutoCompressVCFReader = VCFReader<BufReader<autocompress::Decoder<BufReader<Box<dyn Read>>>>>;

pub fn open_vcf_from_path(path: Option<&str>) -> Result<AutoCompressVCFReader, VCFUtilsError> {
    VCFReader::new(BufReader::new(autocompress::open_or_stdin(path)?)).map_err(|e| e.into())
}

pub fn tsv_reader_builder() -> csv::ReaderBuilder {
    let mut builder = csv::ReaderBuilder::new();
    builder.quoting(false).delimiter(b'\t').escape(None);
    builder
}

pub fn auto_csv_reader_from_path(path: &str, use_header: bool) -> csv::Result<csv::Reader<File>> {
    if path.ends_with(".csv") {
        csv::ReaderBuilder::new()
            .has_headers(use_header)
            .from_path(path)
    } else {
        tsv_reader_builder().has_headers(use_header).from_path(path)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mapping {
    pub mapping: HashMap<U8Vec, U8Vec>,
    pub key_order: Vec<U8Vec>,
    pub value_order: Vec<U8Vec>,
}

pub fn load_mapping<R: Read>(reader: CSVReader<R>) -> Result<Mapping, VCFUtilsError> {
    let mut mapping = HashMap::<U8Vec, U8Vec>::new();
    let mut key_order = Vec::new();
    let mut value_order = Vec::new();

    for row in reader.into_byte_records() {
        let row = row?;
        let key = row
            .get(0)
            .ok_or(VCFUtilsErrorKind::OtherError("no first column"))?;
        let value = row
            .get(1)
            .ok_or(VCFUtilsErrorKind::OtherError("no second column"))?;
        mapping.insert(key.to_vec(), value.to_vec());
        key_order.push(key.to_vec());
        value_order.push(value.to_vec());
    }

    Ok(Mapping {
        mapping,
        key_order,
        value_order,
    })
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_load_mapping() -> Result<(), VCFUtilsError> {
        let mapping = load_mapping(auto_csv_reader_from_path(
            "./testfiles/contig-mapping.csv",
            false,
        )?)?;
        assert_eq!(
            mapping,
            Mapping {
                mapping: [
                    (b"13".to_vec(), b"chr13".to_vec()),
                    (b"14".to_vec(), b"chr14".to_vec())
                ]
                .iter()
                .cloned()
                .collect(),
                key_order: vec![b"13".to_vec(), b"14".to_vec()],
                value_order: vec![b"chr13".to_vec(), b"chr14".to_vec()],
            }
        );
        Ok(())
    }

    #[test]
    fn test_open_vcf_from_path() -> Result<(), VCFUtilsError> {
        open_vcf_from_path(Some("./testfiles/simple1.vcf"))?;
        open_vcf_from_path(Some("./testfiles/1kGP-subset.vcf.gz"))?;
        Ok(())
    }
}
