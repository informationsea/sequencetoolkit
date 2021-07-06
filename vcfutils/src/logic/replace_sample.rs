use crate::error::VCFUtilsError;
use std::collections::HashMap;
use std::hash::BuildHasher;
use std::io::{BufRead, Write};
use vcf::{U8Vec, VCFHeader, VCFReader, VCFRecord, VCFWriter};

pub fn replace_sample<R: BufRead, W: Write, S: BuildHasher>(
    vcf_reader: &mut VCFReader<R>,
    writer: W,
    sample_mapping: &HashMap<U8Vec, U8Vec, S>,
    sample_order: &[U8Vec],
) -> Result<(), VCFUtilsError> {
    let original_sample_index: HashMap<_, _> = vcf_reader
        .header()
        .samples()
        .iter()
        .enumerate()
        .map(|(a, b)| (b.clone(), a))
        .collect();
    let new_sample_index: HashMap<_, _> = sample_order
        .iter()
        .enumerate()
        .map(|(a, b)| (b.clone(), a))
        .collect();
    let mut index_mapping: Vec<_> = sample_mapping
        .iter()
        .filter_map(|(k, v)| {
            if let Some(k_i) = original_sample_index.get(k) {
                if let Some(v_i) = new_sample_index.get(v) {
                    return Some((k_i, v_i, v.clone()));
                }
            }
            None
        })
        .collect();
    index_mapping.sort_by_key(|x| x.1);
    let new_header = VCFHeader::new(
        vcf_reader.header().items().to_vec(),
        index_mapping.iter().map(|x| x.2.clone()).collect(),
    );
    let mut vcf_writer = VCFWriter::new(writer, &new_header)?;
    let mut record = VCFRecord::new(vcf_reader.header().clone());
    while vcf_reader.next_record(&mut record)? {
        let new_genotype: Vec<_> = index_mapping
            .iter()
            .map(|(x, _, _)| record.genotype.get(**x).cloned().unwrap_or_default())
            .collect();
        record.genotype = new_genotype;
        vcf_writer.write_record(&record)?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_replace_sample() -> Result<(), VCFUtilsError> {
        let mapping: HashMap<_, _> = [
            (
                b"ERP001775_HiSeq2000_SAMEA1531955-1".to_vec(),
                b"A".to_vec(),
            ),
            (
                b"ERP001775_HiSeq2000_SAMEA1531955-2".to_vec(),
                b"B".to_vec(),
            ),
        ]
        .iter()
        .cloned()
        .collect();
        let order = vec![b"B".to_vec(), b"A".to_vec()];
        let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
        let mut vcf_reader = VCFReader::new(&vcf_data[..])?;
        let mut write_data: Vec<u8> = Vec::new();
        replace_sample(&mut vcf_reader, &mut write_data, &mapping, &order)?;

        let mut vcf_reader2 = VCFReader::new(&write_data[..])?;
        assert_eq!(
            vcf_reader2.header().samples(),
            &[b"B".to_vec(), b"A".to_vec()]
        );
        let mut record = VCFRecord::new(vcf_reader2.header().clone());
        assert!(vcf_reader2.next_record(&mut record)?);
        assert_eq!(record.genotype(b"A", b"DP"), Some(&vec![b"14".to_vec()]));
        assert_eq!(record.genotype(b"B", b"DP"), Some(&vec![b"19".to_vec()]));

        Ok(())
    }

    #[test]
    fn test_replace_sample2() -> Result<(), VCFUtilsError> {
        let mapping: HashMap<_, _> = [
            (b"SRP150637__HG00099".to_vec(), b"A".to_vec()),
            (b"SRP150637__HG00104".to_vec(), b"B".to_vec()),
            (b"SRP150637__HG00106".to_vec(), b"C".to_vec()),
        ]
        .iter()
        .cloned()
        .collect();
        let order = vec![b"C".to_vec(), b"A".to_vec()];
        let vcf_data = include_bytes!("../../testfiles/1kGP-subset.vcf");
        let mut vcf_reader = VCFReader::new(&vcf_data[..])?;
        let mut write_data: Vec<u8> = Vec::new();
        replace_sample(&mut vcf_reader, &mut write_data, &mapping, &order)?;

        let mut vcf_reader2 = VCFReader::new(&write_data[..])?;
        assert_eq!(
            vcf_reader2.header().samples(),
            &[b"C".to_vec(), b"A".to_vec()]
        );
        let mut record = VCFRecord::new(vcf_reader2.header().clone());
        assert!(vcf_reader2.next_record(&mut record)?);
        assert_eq!(record.genotype(b"A", b"DP"), Some(&vec![b"31".to_vec()]));
        assert_eq!(record.genotype(b"C", b"DP"), Some(&vec![b"26".to_vec()]));

        Ok(())
    }
}
