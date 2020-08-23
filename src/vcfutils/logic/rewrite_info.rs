use crate::vcfutils::error::VCFUtilsError;
use std::collections::HashMap;
use std::io;

pub fn rewrite_info<R: io::BufRead, W: io::Write>(
    reader: &mut vcf::VCFReader<R>,
    writer: W,
    info_list: &[vcf::U8Vec],
) -> Result<(), VCFUtilsError> {
    let mut new_header_items: Vec<_> = reader
        .header()
        .items()
        .iter()
        .filter(|x| match x.contents() {
            vcf::VCFHeaderContent::INFO { .. } => false,
            _ => true,
        })
        .cloned()
        .collect();
    let info_line_map: HashMap<_, _> = reader
        .header()
        .items()
        .iter()
        .filter_map(|x| match x.contents() {
            vcf::VCFHeaderContent::INFO { id, .. } => {
                if info_list.contains(id) {
                    Some((id.clone(), x.clone()))
                } else {
                    None
                }
            }
            _ => None,
        })
        .collect();
    new_header_items.append(
        &mut info_list
            .iter()
            .filter_map(|x| info_line_map.get(x))
            .cloned()
            .collect(),
    );
    let new_header = vcf::VCFHeader::new(new_header_items, reader.header().samples().to_vec());
    let mut vcf_writer = vcf::VCFWriter::new(writer, &new_header)?;
    let mut record = vcf::VCFRecord::new(new_header);

    while reader.next_record(&mut record)? {
        let new_info: Vec<_> = info_list
            .iter()
            .filter_map(|x| record.info(x).map(|y| (x.clone(), y.clone())))
            .collect();
        record.info = new_info;
        vcf_writer.write_record(&record)?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::HashSet;
    use std::io;
    #[test]
    fn test_rewrite_info() -> Result<(), VCFUtilsError> {
        let mut write_result: Vec<u8> = Vec::new();
        let mut vcf_reader = vcf::VCFReader::new(io::BufReader::new(
            &include_bytes!("../../../testfiles/vcfutils/simple1.vcf")[..],
        ))?;

        rewrite_info(
            &mut vcf_reader,
            &mut write_result,
            &[b"AF".to_vec(), b"AN".to_vec()][..],
        )?;

        let mut vcf_reader2 = vcf::VCFReader::new(io::BufReader::new(&write_result[..]))?;
        assert_eq!(
            [b"AF".to_vec(), b"AN".to_vec()]
                .iter()
                .cloned()
                .collect::<HashSet<_>>(),
            vcf_reader2
                .header()
                .info_list()
                .map(|x| x.to_vec())
                .collect::<HashSet<_>>()
        );
        let mut record = vcf::VCFRecord::new(vcf_reader2.header().clone());
        while vcf_reader2.next_record(&mut record)? {
            assert_eq!(
                record.info.iter().map(|x| &x.0).collect::<Vec<_>>(),
                vec![&b"AF".to_vec(), &b"AN".to_vec()]
            );
        }
        Ok(())
    }
}
