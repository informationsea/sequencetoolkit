use crate::vcfutils::error::VCFUtilsError;
use std::collections::{HashMap, HashSet};
use std::io;

pub fn rewrite_format<R: io::BufRead, W: io::Write, S: std::hash::BuildHasher>(
    reader: &mut vcf::VCFReader<R>,
    writer: W,
    format_list: &HashSet<vcf::U8Vec, S>,
) -> Result<(), VCFUtilsError> {
    let mut new_header_items: Vec<_> = reader
        .header()
        .items()
        .iter()
        .filter(|x| match x.contents() {
            vcf::VCFHeaderContent::FORMAT { .. } => false,
            _ => true,
        })
        .cloned()
        .collect();
    let info_line_map: HashMap<_, _> = reader
        .header()
        .items()
        .iter()
        .filter_map(|x| match x.contents() {
            vcf::VCFHeaderContent::FORMAT { id, .. } => {
                if format_list.contains(id) {
                    Some((id.clone(), x.clone()))
                } else {
                    None
                }
            }
            _ => None,
        })
        .collect();
    new_header_items.append(
        &mut format_list
            .iter()
            .filter_map(|x| info_line_map.get(x))
            .cloned()
            .collect(),
    );
    let new_header = vcf::VCFHeader::new(new_header_items, reader.header().samples().to_vec());
    let mut vcf_writer = vcf::VCFWriter::new(writer, &new_header)?;
    let mut record = vcf::VCFRecord::new(new_header);

    while reader.next_record(&mut record)? {
        let new_format: Vec<_> = record
            .format
            .iter()
            .map(|x| format_list.contains(x))
            .collect();
        record.format = record
            .format
            .into_iter()
            .zip(new_format.iter())
            .filter(|(_, b)| **b)
            .map(|(x, _)| x)
            .collect();
        for one in record.genotype.iter_mut() {
            *one = (one)
                .iter()
                .zip(new_format.iter())
                .filter(|(_, b)| **b)
                .map(|(x, _)| x.clone())
                .collect();
        }

        vcf_writer.write_record(&record)?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use std::io;
    #[test]
    fn test_rewrite_format() -> Result<(), VCFUtilsError> {
        let mut write_result: Vec<u8> = Vec::new();
        let mut vcf_reader = vcf::VCFReader::new(io::BufReader::new(
            &include_bytes!("../../../testfiles/vcfutils/simple1.vcf")[..],
        ))?;

        rewrite_format::<_, _, std::collections::hash_map::RandomState>(
            &mut vcf_reader,
            &mut write_result,
            &[b"GT".to_vec(), b"DP".to_vec()].iter().cloned().collect(),
        )?;

        let mut vcf_reader2 = vcf::VCFReader::new(io::BufReader::new(&write_result[..]))?;
        assert_eq!(
            [b"GT".to_vec(), b"DP".to_vec()]
                .iter()
                .cloned()
                .collect::<HashSet<_>>(),
            vcf_reader2
                .header()
                .format_list()
                .map(|x| x.to_vec())
                .collect::<HashSet<_>>()
        );
        let mut record = vcf::VCFRecord::new(vcf_reader2.header().clone());
        while vcf_reader2.next_record(&mut record)? {
            assert_eq!(
                record.format.iter().collect::<Vec<_>>(),
                vec![&b"GT".to_vec(), &b"DP".to_vec()]
            );
            for one_genotype in &record.genotype {
                assert_eq!(one_genotype.len(), 2);
                assert_eq!(one_genotype[0][0].len(), 3);
            }
        }

        Ok(())
    }
}
