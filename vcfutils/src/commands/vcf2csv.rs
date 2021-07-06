use super::Command;
use crate::error::VCFUtilsError;
use crate::logic::vcf2table::{
    create_header_line, vcf2table, vcf2table_set_data_type, VCF2CSVConfig,
};
use crate::utils;
use crate::utils::tablewriter::{CSVWriter, TSVWriter, TableWriter, XlsxSheetWriter};
use clap::{App, Arg, ArgMatches};
use std::boxed::Box;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use vcf::{VCFHeader, VCFReader};

pub struct VCF2CSV;

impl Command for VCF2CSV {
    fn command_name(&self) -> &'static str {
        "vcf2csv"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Convert VCF to csv")
            .arg(
                Arg::with_name("input")
                    .index(1)
                    .takes_value(true)
                    .help("Input VCF file"),
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("output")
                    .takes_value(true)
                    .help("Output file"),
            )
            .arg(
                Arg::with_name("canonical-list")
                    .short("c")
                    .long("canonical-list")
                    .takes_value(true)
                    .help("Canonical transcript list (created with extract-canonical command)"),
            )
            .arg(
                Arg::with_name("datatype")
                    .short("t")
                    .long("data-type")
                    .alias("datatype")
                    .takes_value(true)
                    .help("Data type to output")
                    .possible_values(&["auto", "xlsx", "csv", "tsv"])
                    .default_value("auto"),
            )
            .arg(
                Arg::with_name("split-multi-allelic")
                    .short("m")
                    .long("split-multi-allelic")
                    .help("Split multi allelic sites into multiple lines"),
            )
            .arg(
                Arg::with_name("decode-genotype")
                    .short("d")
                    .long("decode-genotype")
                    .help("Decode GT format tag number into alleles"),
            )
            .arg(
                Arg::with_name("info")
                    .short("i")
                    .long("info")
                    .help("INFO tags to include")
                    .takes_value(true)
                    .multiple(true),
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("FORMAT tags to include")
                    .takes_value(true)
                    .multiple(true),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        let output_type = match matches.value_of("datatype").unwrap() {
            "csv" => "csv",
            "tsv" => "tsv",
            "xlsx" => "xlsx",
            "auto" => {
                if let Some(output_path) = matches.value_of("output") {
                    if output_path.ends_with(".xlsx") {
                        "xlsx"
                    } else if output_path.ends_with(".csv") {
                        "csv"
                    } else {
                        "tsv"
                    }
                } else {
                    "tsv"
                }
            }
            _ => unreachable!(),
        };

        let mut vcf_reader = utils::open_vcf_from_path(matches.value_of("input"))?;
        let mut writer: Box<dyn TableWriter> = match output_type {
            "csv" => Box::new(CSVWriter::new(autocompress::create_or_stdout(
                matches.value_of("output"),
            )?)),
            "tsv" => Box::new(TSVWriter::new(autocompress::create_or_stdout(
                matches.value_of("output"),
            )?)),
            "xlsx" => {
                return Ok(run_xlsx_mode(vcf_reader, matches)?);
            }
            _ => unreachable!(),
        };
        let config = create_config(&vcf_reader.header(), matches)?;
        let header_contents = create_header_line(&vcf_reader.header(), &config);
        vcf2table(&mut vcf_reader, &header_contents, &config, &mut *writer)?;
        Ok(())
    }
}

pub fn create_config(
    header: &VCFHeader,
    matches: &ArgMatches<'static>,
) -> Result<VCF2CSVConfig, VCFUtilsError> {
    let canonical_list = if let Some(canonical_file) = matches.value_of("canonical-list") {
        let mut canonical_list = HashSet::new();

        let mut reader = BufReader::new(autocompress::open(canonical_file)?);
        let mut buffer = Vec::new();

        while reader.read_until(b'\n', &mut buffer)? > 0 {
            if buffer.ends_with(b"\r\n") {
                buffer.pop();
                buffer.pop();
            } else {
                buffer.pop();
            }
            canonical_list.insert(buffer.to_vec());
            buffer.clear();
        }

        Some(canonical_list)
    } else {
        None
    };

    Ok(VCF2CSVConfig {
        split_multi_allelic: matches.is_present("split-multi-allelic"),
        decoded_genotype: matches.is_present("decode-genotype"),
        canonical_list,
        info_list: matches
            .values_of("info")
            .map(|x| {
                x.map(|y| y.as_bytes().to_vec())
                    .collect::<Vec<vcf::U8Vec>>()
            })
            .unwrap_or_else(|| {
                header
                    .items()
                    .iter()
                    .filter_map(|x| match x.contents() {
                        vcf::VCFHeaderContent::INFO { id, .. } => Some(id.to_vec()),
                        _ => None,
                    })
                    .collect()
            }),
        format_list: matches
            .values_of("format")
            .map(|x| {
                x.map(|y| y.as_bytes().to_vec())
                    .collect::<Vec<vcf::U8Vec>>()
            })
            .unwrap_or_else(|| {
                header
                    .items()
                    .iter()
                    .filter_map(|x| match x.contents() {
                        vcf::VCFHeaderContent::FORMAT { id, .. } => Some(id.to_vec()),
                        _ => None,
                    })
                    .collect()
            }),
    })
}

fn run_xlsx_mode<R: BufRead>(
    mut vcf_reader: VCFReader<R>,
    matches: &ArgMatches<'static>,
) -> Result<(), VCFUtilsError> {
    let workbook = xlsxwriter::Workbook::new(
        matches
            .value_of("output")
            .expect("Output path is required for xlsx output mode"),
    );
    let mut sheet = workbook.add_worksheet(None)?;
    let mut writer = XlsxSheetWriter::new(&mut sheet);
    let config = create_config(&vcf_reader.header(), matches)?;
    let header_contents = create_header_line(&vcf_reader.header(), &config);
    vcf2table_set_data_type(&header_contents, &mut writer)?;
    let row = vcf2table(&mut vcf_reader, &header_contents, &config, &mut writer)?;
    sheet.autofilter(0, 0, row, (header_contents.len() - 1) as u16)?;
    sheet.freeze_panes(1, 0);
    workbook.close()?;
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::error::VCFUtilsError;
    #[test]
    fn test_create_config() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
        let vcf_reader = VCFReader::new(&vcf_data[..])?;
        let matches = VCF2CSV
            .config_subcommand(App::new("test"))
            .get_matches_from_safe(["prog", "-t", "xlsx", "-c", "testfiles/canonical.txt"].iter())
            .unwrap();
        let config = create_config(&vcf_reader.header(), &matches)?;
        assert_eq!(config.split_multi_allelic, false);
        assert_eq!(config.decoded_genotype, false);
        assert_eq!(
            config.info_list,
            vec![
                b"AC".to_vec(),
                b"AF".to_vec(),
                b"AN".to_vec(),
                b"DP".to_vec()
            ]
        );
        assert_eq!(
            config.format_list,
            vec![b"AD".to_vec(), b"DP".to_vec(), b"GT".to_vec()]
        );
        assert!(config
            .canonical_list
            .as_ref()
            .unwrap()
            .contains(&b"ENST00000384061.1".to_vec()));
        assert!(config
            .canonical_list
            .as_ref()
            .unwrap()
            .contains(&b"ENST00000533876.1_1".to_vec()));

        let matches = VCF2CSV
            .config_subcommand(App::new("test"))
            .get_matches_from_safe(
                ["prog", "-m", "-i", "AC", "-f", "AD", "--decode-genotype"].iter(),
            )
            .unwrap();
        let config = create_config(&vcf_reader.header(), &matches)?;
        assert_eq!(
            config,
            VCF2CSVConfig {
                split_multi_allelic: true,
                decoded_genotype: true,
                canonical_list: None,
                info_list: vec![b"AC".to_vec(),],
                format_list: vec![b"AD".to_vec()]
            }
        );
        Ok(())
    }
}
