use crate::vcfutils::error::VCFUtilsError;
use std::collections::{HashMap, HashSet};
use std::io::prelude::*;
use vcf::U8Vec;

pub fn scan_contig(reader: &mut impl BufRead) -> Result<HashSet<U8Vec>, VCFUtilsError> {
    let mut contigs = HashSet::new();
    let mut buffer = Vec::new();

    loop {
        buffer.clear();
        let read_size = reader.read_until(b'\n', &mut buffer)?;
        if read_size == 0 {
            break;
        }
        if buffer[0] == b'#' {
            continue;
        }
        if let Some(x) = buffer.split(|x| *x == b'\t').next() {
            if !contigs.contains(x) {
                contigs.insert(x.to_vec());
            }
        }
    }

    Ok(contigs)
}

pub fn add_contig(
    reader: &mut impl BufRead,
    writer: &mut impl Write,
    contig: &HashSet<U8Vec>,
) -> Result<(), VCFUtilsError> {
    let mut buffer = Vec::new();
    let mut line = 0;
    let mut contig_header: HashMap<U8Vec, vcf::VCFHeaderLine> = HashMap::new();

    loop {
        line += 1;
        buffer.clear();
        let read_size = reader.read_until(b'\n', &mut buffer)?;
        if read_size == 0 {
            break;
        }
        if buffer.starts_with(b"##") {
            let parsed_header = vcf::VCFHeaderLine::from_bytes(&buffer, line)?;
            match parsed_header.contents() {
                vcf::VCFHeaderContent::Contig { id, .. } => {
                    contig_header.insert(id.to_vec(), parsed_header);
                }
                _ => {
                    writer.write_all(&buffer)?;
                }
            }
        } else {
            let mut contig_list: Vec<_> = contig.iter().collect();
            contig_list.sort();
            for one in contig_list {
                if let Some(item) = contig_header.get(one) {
                    writer.write_all(item.line())?;
                } else {
                    writer.write_all(b"##contig=<ID=")?;
                    writer.write_all(one)?;
                    writer.write_all(b">\n")?;
                }
            }
            writer.write_all(&buffer)?;
            break;
        }
    }

    std::io::copy(reader, writer)?;

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_add_contig() -> Result<(), VCFUtilsError> {
        let sample_vcf = br#"##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##source=CombineGVCFs
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SRP150637__HG00099	SRP150637__HG00100	SRP150637__HG00102	SRP150637__HG00104	SRP150637__HG00106
13	32872836	.	A	C	495.23	.	AC=1;AF=0.1;AN=10;BaseQRankSum=-1.425;DP=150;ExcessHet=3.0103;FS=5.469;MLEAC=1;MLEAF=0.1;MQ=60;MQRankSum=0;QD=16.51;ReadPosRankSum=0.381;SOR=0.589	GT:AD:DP:GQ:PL	0/0:31,0:31:82:0,82,1043	0/1:12,18:30:99:505,0,327	0/0:28,0:28:84:0,84,933	0/0:35,0:35:99:0,99,1485	0/0:26,0:26:70:0,70,864
"#;

        let contigs: HashSet<_> = [b"1".to_vec(), b"10".to_vec()].iter().cloned().collect();
        let mut buffer = Vec::<u8>::new();

        add_contig(&mut BufReader::new(&sample_vcf[..]), &mut buffer, &contigs)?;

        let expected_vcf = br#"##fileformat=VCFv4.2
##source=CombineGVCFs
##contig=<ID=1,length=249250621>
##contig=<ID=10>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SRP150637__HG00099	SRP150637__HG00100	SRP150637__HG00102	SRP150637__HG00104	SRP150637__HG00106
13	32872836	.	A	C	495.23	.	AC=1;AF=0.1;AN=10;BaseQRankSum=-1.425;DP=150;ExcessHet=3.0103;FS=5.469;MLEAC=1;MLEAF=0.1;MQ=60;MQRankSum=0;QD=16.51;ReadPosRankSum=0.381;SOR=0.589	GT:AD:DP:GQ:PL	0/0:31,0:31:82:0,82,1043	0/1:12,18:30:99:505,0,327	0/0:28,0:28:84:0,84,933	0/0:35,0:35:99:0,99,1485	0/0:26,0:26:70:0,70,864
"#;

        assert_eq!(&expected_vcf[..], &buffer[..]);

        Ok(())
    }
    #[test]
    fn test_scan_contig() -> Result<(), VCFUtilsError> {
        let mut reader = BufReader::new(autocompress::open(
            "testfiles/vcfutils/SRP225680__ISMMS_pt_91_WES_TUMOR.hs37d5-filtered.snpeff.vcf.gz",
        )?);
        let contigs = scan_contig(&mut reader)?;

        let expected: HashSet<_> = [
            b"1".to_vec(),
            b"10".to_vec(),
            b"11".to_vec(),
            b"12".to_vec(),
            b"13".to_vec(),
            b"14".to_vec(),
            b"15".to_vec(),
            b"16".to_vec(),
            b"17".to_vec(),
            b"18".to_vec(),
            b"19".to_vec(),
            b"2".to_vec(),
            b"20".to_vec(),
            b"21".to_vec(),
            b"22".to_vec(),
            b"3".to_vec(),
            b"4".to_vec(),
            b"5".to_vec(),
            b"6".to_vec(),
            b"7".to_vec(),
            b"8".to_vec(),
            b"9".to_vec(),
            b"X".to_vec(),
        ]
        .iter()
        .cloned()
        .collect();

        assert_eq!(contigs, expected);

        Ok(())
    }
}
