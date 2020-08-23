use crate::vcfutils::error::VCFUtilsError;
use std::collections::HashMap;
use std::hash::BuildHasher;
use std::io::prelude::*;

pub fn replace_contig<W: Write, R: BufRead, S: BuildHasher>(
    mut reader: R,
    mut writer: W,
    mapping: &HashMap<vcf::U8Vec, vcf::U8Vec, S>,
) -> Result<(), VCFUtilsError> {
    let mut line = Vec::new();
    let mut count = 0;
    while reader.read_until(b'\n', &mut line)? > 0 {
        count += 1;
        if line.starts_with(b"##contig") {
            let parsed_header = vcf::VCFHeaderLine::from_bytes(&line, count)?;
            if let vcf::VCFHeaderContent::Contig { id, length } = parsed_header.contents() {
                writer.write_all(b"##contig=<ID=")?;
                writer.write_all(mapping.get(id).unwrap_or(id))?;
                if let Some(length) = length {
                    writer.write_all(b",length=")?;
                    write!(writer, "{}", length)?;
                }
                writer.write_all(b">\n")?;
            } else {
                unreachable!()
            }
        } else if line.starts_with(b"#") {
            writer.write_all(&line)?;
        } else {
            for (i, data) in line.splitn(2, |x| *x == b'\t').enumerate() {
                if i == 0 {
                    writer.write_all(mapping.get(data).map(|x| -> &[u8] { x }).unwrap_or(data))?;
                    writer.write_all(b"\t")?;
                } else {
                    writer.write_all(data)?;
                }
            }
        }
        line.clear();
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::HashSet;
    use vcf::{U8Vec, VCFHeaderContent, VCFReader, VCFRecord};

    #[test]
    fn test_replace_contig() -> Result<(), VCFUtilsError> {
        let mapping: HashMap<U8Vec, U8Vec> = [(b"13".to_vec(), b"chr13".to_vec())]
            .iter()
            .cloned()
            .collect();
        let input_vcf = include_bytes!("../../../testfiles/vcfutils/simple1.vcf");
        let mut result = Vec::new();
        replace_contig(&input_vcf[..], &mut result, &mapping)?;

        std::fs::File::create("target/replace.vcf")?.write_all(&result)?;

        let mut vcf_reader = VCFReader::new(&result[..])?;
        let contig: HashSet<U8Vec> = vcf_reader
            .header()
            .items()
            .iter()
            .filter_map(|x| match x.contents() {
                VCFHeaderContent::Contig { id, .. } => Some(id.clone()),
                _ => None,
            })
            .collect();
        assert_eq!(
            contig,
            [b"chr13".to_vec(), b"14".to_vec()]
                .iter()
                .cloned()
                .collect::<HashSet<_>>()
        );
        let mut record = VCFRecord::new(vcf_reader.header().clone());
        assert!(vcf_reader.next_record(&mut record)?);
        assert_eq!(record.chromosome, b"chr13");
        assert!(vcf_reader.next_record(&mut record)?);
        assert_eq!(record.chromosome, b"chr13");
        assert!(vcf_reader.next_record(&mut record)?);
        assert_eq!(record.chromosome, b"chr13");

        Ok(())
    }
}
