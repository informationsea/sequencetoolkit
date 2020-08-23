use super::vcf2table::{create_header_line, HeaderType, VCF2CSVConfig};
use crate::vcfutils::error::VCFUtilsError;
use itertools::join;
use nom::character::is_alphanumeric;
use std::string::String;
use vcf::{self, VCFHeader, ValueType};

pub fn generate_sql(
    header: &VCFHeader,
    config: &VCF2CSVConfig,
    table_name: &str,
) -> Result<String, VCFUtilsError> {
    let header_items = create_header_line(header, config);
    let mut sql = format!("CREATE TABLE {} (\n", table_name);

    sql.push_str(&join(
        header_items
            .iter()
            .filter(|x| **x != HeaderType::AltIndex && **x != HeaderType::VcfLine)
            .map(|x| {
                format!(
                    "    {} {}",
                    String::from_utf8(clean_name(x.to_string().as_bytes())).unwrap(),
                    match x {
                        HeaderType::VcfLine => "",
                        HeaderType::AltIndex => "",
                        HeaderType::GeneChange => "TEXT",
                        HeaderType::GeneName => "TEXT",
                        HeaderType::TranscriptName => "TEXT",
                        HeaderType::AminoChange => "TEXT",
                        HeaderType::CDSChange => "TEXT",
                        HeaderType::CHROM => "TEXT",
                        HeaderType::POS => "INTEGER",
                        HeaderType::ID => "TEXT",
                        HeaderType::REF => "TEXT",
                        HeaderType::ALT => "TEXT",
                        HeaderType::QUAL => "REAL",
                        HeaderType::FILTER => "TEXT",
                        HeaderType::Info(_, _, value_type, ..)
                        | HeaderType::Genotype(_, _, _, value_type, ..) => match value_type {
                            ValueType::Integer => "INTEGER",
                            ValueType::Float => "REAL",
                            _ => "TEXT",
                        },
                        HeaderType::SnpEffImpact(_) => "TEXT",
                        HeaderType::SnpEff => "TEXT",
                    }
                )
            }),
        ",\n",
    ));
    sql.push_str("\n);");
    Ok(sql)
}

fn clean_name(name: &[u8]) -> Vec<u8> {
    name.iter()
        .map(|x| if is_alphanumeric(*x) { *x } else { b'_' })
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_clean_name() {
        assert_eq!(clean_name(b"abc123 !-=_"), b"abc123_____");
    }

    #[test]
    fn test_generate_sql() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../../testfiles/vcfutils/simple1.vcf");
        let config = VCF2CSVConfig {
            split_multi_allelic: true,
            decoded_genotype: false,
            canonical_list: None,
            info_list: vec![b"AC".to_vec(), b"AN".to_vec(), b"AF".to_vec()],
            format_list: vec![b"GT".to_vec(), b"AD".to_vec(), b"DP".to_vec()],
        };
        let vcf_reader = vcf::VCFReader::new(&vcf_data[..])?;
        assert_eq!(
            generate_sql(&vcf_reader.header(), &config, "vcf_table")?,
            r#"CREATE TABLE vcf_table (
    CHROM TEXT,
    POS INTEGER,
    ID TEXT,
    REF TEXT,
    ALT TEXT,
    QUAL REAL,
    FILTER TEXT,
    AC INTEGER,
    AN INTEGER,
    AF REAL,
    ERP001775_HiSeq2000_SAMEA1531955_1__GT TEXT,
    ERP001775_HiSeq2000_SAMEA1531955_1__AD__Ref INTEGER,
    ERP001775_HiSeq2000_SAMEA1531955_1__AD__Alt INTEGER,
    ERP001775_HiSeq2000_SAMEA1531955_1__DP INTEGER,
    ERP001775_HiSeq2000_SAMEA1531955_2__GT TEXT,
    ERP001775_HiSeq2000_SAMEA1531955_2__AD__Ref INTEGER,
    ERP001775_HiSeq2000_SAMEA1531955_2__AD__Alt INTEGER,
    ERP001775_HiSeq2000_SAMEA1531955_2__DP INTEGER
);"#
        );

        Ok(())
    }
}
