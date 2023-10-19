use crate::error::VCFUtilsError;
use crate::logic::vcf2table::{
    create_header_line, merge_header_contents, vcf2table, vcf2table_set_data_type, VCF2CSVConfig,
};
use crate::utils;
use crate::utils::tablewriter::{CSVWriter, TSVWriter, TableWriter, XlsxSheetWriter};
use clap::{Args, ValueEnum};
use std::boxed::Box;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use vcf::VCFHeader;

pub trait TableConfig {
    fn canonical_list(&self) -> Option<&str>;
    fn input(&self) -> &[String];
    fn group_names(&self) -> Option<&[String]>;
    fn split_multi_allelic(&self) -> bool;
    fn decode_genotype(&self) -> bool;
    fn info_list(&self) -> Option<&[String]>;
    fn format_list(&self) -> Option<&[String]>;
    fn replace_sample_name(&self) -> Option<&[String]>;
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum DataType {
    Auto,
    CSV,
    TSV,
    XLSX,
}

#[derive(Args, Debug)]
#[command(about = "Convert VCF to csv", version, author)]
pub struct VCF2CSV {
    #[arg(help = "Input VCF files")]
    input: Vec<String>,
    #[arg(short, long, help = "Output file")]
    output: Option<String>,
    #[arg(
        short,
        long,
        help = "Canonical transcript list (created with extract-canonical command)"
    )]
    canonical_list: Option<String>,
    #[arg(
        short = 't',
        long = "type",
        help = "Data type to output",
        default_value = "auto"
    )]
    datatype: DataType,
    #[arg(
        short = 'm',
        long,
        help = "Split multi allelic sites into multiple lines"
    )]
    split_multi_allelic: bool,
    #[arg(short = 'd', long, help = "Decode GT format tag number into alleles")]
    decode_genotype: bool,
    #[arg(short = 'i', long, help = "INFO tags to include")]
    info: Option<Vec<String>>,
    #[arg(short = 'f', long, help = "FORMAT tags to include")]
    format: Option<Vec<String>>,
    #[arg(
        short = 'g',
        long,
        help = "List of group names for each input VCF file."
    )]
    group_names: Option<Vec<String>>,
    #[arg(short = 'r', long, help = "Replace sample name for each input vcf")]
    replace_sample_name: Option<Vec<String>>,
}

impl TableConfig for VCF2CSV {
    fn canonical_list(&self) -> Option<&str> {
        self.canonical_list.as_deref()
    }
    fn input(&self) -> &[String] {
        self.input.as_slice()
    }
    fn replace_sample_name(&self) -> Option<&[String]> {
        self.replace_sample_name.as_deref()
    }
    fn group_names(&self) -> Option<&[String]> {
        self.group_names.as_deref()
    }
    fn split_multi_allelic(&self) -> bool {
        self.split_multi_allelic
    }
    fn decode_genotype(&self) -> bool {
        self.decode_genotype
    }
    fn info_list(&self) -> Option<&[String]> {
        self.info.as_deref()
    }
    fn format_list(&self) -> Option<&[String]> {
        self.format.as_deref()
    }
}

impl VCF2CSV {
    // fn command_name(&self) -> &'static str {
    //     "vcf2csv"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("Convert VCF to csv")
    //         .arg(
    //             Arg::with_name("input")
    //                 .index(1)
    //                 .takes_value(true)
    //                 .multiple(true)
    //                 .help("Input VCF file"),
    //         )
    //         .arg(
    //             Arg::with_name("output")
    //                 .short("o")
    //                 .long("output")
    //                 .takes_value(true)
    //                 .help("Output file"),
    //         )
    //         .arg(
    //             Arg::with_name("canonical-list")
    //                 .short("c")
    //                 .long("canonical-list")
    //                 .takes_value(true)
    //                 .help("Canonical transcript list (created with extract-canonical command)"),
    //         )
    //         .arg(
    //             Arg::with_name("datatype")
    //                 .short("t")
    //                 .long("data-type")
    //                 .alias("datatype")
    //                 .takes_value(true)
    //                 .help("Data type to output")
    //                 .possible_values(&["auto", "xlsx", "csv", "tsv"])
    //                 .default_value("auto"),
    //         )
    //         .arg(
    //             Arg::with_name("split-multi-allelic")
    //                 .short("m")
    //                 .long("split-multi-allelic")
    //                 .help("Split multi allelic sites into multiple lines"),
    //         )
    //         .arg(
    //             Arg::with_name("decode-genotype")
    //                 .short("d")
    //                 .long("decode-genotype")
    //                 .help("Decode GT format tag number into alleles"),
    //         )
    //         .arg(
    //             Arg::with_name("info")
    //                 .short("i")
    //                 .long("info")
    //                 .help("INFO tags to include")
    //                 .takes_value(true)
    //                 .multiple(true),
    //         )
    //         .arg(
    //             Arg::with_name("format")
    //                 .short("f")
    //                 .long("format")
    //                 .help("FORMAT tags to include")
    //                 .takes_value(true)
    //                 .multiple(true),
    //         )
    //         .arg(
    //             Arg::with_name("group-names")
    //                 .short("g")
    //                 .long("group-names")
    //                 .takes_value(true)
    //                 .multiple(true)
    //                 .help("add Group Name column and fill with a value for each input vcf"),
    //         )
    //         .arg(
    //             Arg::with_name("replace-sample-name")
    //                 .short("r")
    //                 .help("Replace sample name for each input vcf")
    //                 .takes_value(true)
    //                 .multiple(true),
    //         )
    // }

    pub fn run(&self) -> anyhow::Result<()> {
        let output_type = match self.datatype {
            DataType::CSV => "csv",
            DataType::TSV => "tsv",
            DataType::XLSX => "xlsx",
            DataType::Auto => {
                if let Some(output_path) = self.output.as_ref() {
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
        };

        // let vcf_files: Vec<_> = self
        //     .input
        //     .iter()
        //     .map(|x| x.collect())
        //     .unwrap_or_else(|| vec![]);

        if output_type == "xlsx" {
            return Ok(self.run_xlsx_mode(&self.input)?);
        }

        let mut vcf_reader = utils::open_vcf_from_path(self.input.get(0).map(|x| x.as_str()))?;
        let mut writer: Box<dyn TableWriter> = match output_type {
            "csv" => Box::new(CSVWriter::new(autocompress::create_or_stdout(
                self.output.as_deref(),
                autocompress::CompressionLevel::Default,
            )?)),
            "tsv" => Box::new(TSVWriter::new(autocompress::create_or_stdout(
                self.output.as_deref(),
                autocompress::CompressionLevel::Default,
            )?)),
            _ => unreachable!(),
        };
        let config = create_config(&vcf_reader.header(), self)?;
        let header_contents = create_header_line(&vcf_reader.header(), &config);

        vcf2table(
            &mut vcf_reader,
            &header_contents,
            &config,
            config.group_names.as_ref().map(|x| x.get(0)).flatten(),
            true,
            &mut *writer,
        )?;

        if self.input.len() > 1 {
            for (i, one_vcf_name) in self.input[1..].iter().enumerate() {
                let mut vcf_reader = utils::open_vcf_from_path(Some(one_vcf_name))?;
                let config = create_config(&vcf_reader.header(), self)?;
                let new_header_contents = create_header_line(&vcf_reader.header(), &config);
                let merged_header_contents =
                    merge_header_contents(&header_contents, &new_header_contents);
                vcf2table(
                    &mut vcf_reader,
                    &merged_header_contents,
                    &config,
                    config.group_names.as_ref().map(|x| x.get(i + 1)).flatten(),
                    false,
                    &mut writer,
                )?;
            }
        }
        Ok(())
    }

    fn run_xlsx_mode(&self, vcf_inputs: &[String]) -> Result<(), VCFUtilsError> {
        let workbook = xlsxwriter::Workbook::new(
            self.output
                .as_deref()
                .expect("Output path is required for xlsx output mode"),
        )?;

        let mut first_vcf_reader =
            utils::open_vcf_from_path(vcf_inputs.get(0).map(|x| x.as_str()))?;

        let mut sheet = workbook.add_worksheet(None)?;
        let mut writer = XlsxSheetWriter::new(&mut sheet);
        let config = create_config(&first_vcf_reader.header(), self)?;
        let header_contents = create_header_line(&first_vcf_reader.header(), &config);
        vcf2table_set_data_type(&header_contents, &mut writer)?;

        let mut row = vcf2table(
            &mut first_vcf_reader,
            &header_contents,
            &config,
            config.group_names.as_ref().map(|x| x.get(0)).flatten(),
            true,
            &mut writer,
        )?;

        if vcf_inputs.len() > 1 {
            for (i, one_vcf_name) in vcf_inputs[1..].iter().enumerate() {
                let mut vcf_reader = utils::open_vcf_from_path(Some(one_vcf_name))?;
                let config = create_config(&vcf_reader.header(), self)?;
                let new_header_contents = create_header_line(&vcf_reader.header(), &config);
                let merged_header_contents =
                    merge_header_contents(&header_contents, &new_header_contents);
                vcf2table_set_data_type(&merged_header_contents, &mut writer)?;
                let additional_row = vcf2table(
                    &mut vcf_reader,
                    &merged_header_contents,
                    &config,
                    config.group_names.as_ref().map(|x| x.get(i + 1)).flatten(),
                    false,
                    &mut writer,
                )?;
                row += additional_row;
            }
        }
        sheet.autofilter(0, 0, row, (header_contents.len() - 1) as u16)?;
        sheet.freeze_panes(1, 0);
        workbook.close()?;
        Ok(())
    }
}

pub fn create_config(
    header: &VCFHeader,
    matches: &impl TableConfig,
) -> Result<VCF2CSVConfig, VCFUtilsError> {
    let canonical_list = if let Some(canonical_file) = matches.canonical_list() {
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

    if matches.input().len() > 1 {
        if matches.replace_sample_name().is_none() {
            log::error!("--replace-sample-name option is required for multi VCF inputs");
            return Err(VCFUtilsError::OtherError(
                "--replace-sample-name option is required for multi VCF inputs",
            ));
        }
    }

    let group_names = matches
        .group_names()
        .map(|x| x.iter().map(|y| y.as_bytes().to_vec()).collect())
        .or_else(|| {
            if matches.input().len() > 1 {
                Some(
                    matches
                        .input()
                        .iter()
                        .map(|y| y.as_bytes().to_vec())
                        .collect(),
                )
            } else {
                None
            }
        });

    Ok(VCF2CSVConfig {
        split_multi_allelic: matches.split_multi_allelic(),
        decoded_genotype: matches.decode_genotype(),
        canonical_list,
        info_list: matches
            .info_list()
            .map(|x| {
                x.iter()
                    .map(|y| y.as_bytes().to_vec())
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
            .format_list()
            .map(|x| {
                x.iter()
                    .map(|y| y.as_bytes().to_vec())
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
        replace_sample_name: matches
            .replace_sample_name()
            .map(|x| x.iter().map(|y| y.as_bytes().to_vec()).collect()),
        group_names,
    })
}

// #[cfg(test)]
// mod test {
//     use super::*;
//     use crate::error::VCFUtilsError;
//     use vcf::VCFReader;
//     #[test]
//     fn test_create_config() -> Result<(), VCFUtilsError> {
//         let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
//         let vcf_reader = VCFReader::new(&vcf_data[..])?;
//         let matches = VCF2CSV
//             .config_subcommand(App::new("test"))
//             .get_matches_from_safe(["prog", "-t", "xlsx", "-c", "testfiles/canonical.txt"].iter())
//             .unwrap();
//         let config = create_config(&vcf_reader.header(), &matches)?;
//         assert_eq!(config.split_multi_allelic, false);
//         assert_eq!(config.decoded_genotype, false);
//         assert_eq!(
//             config.info_list,
//             vec![
//                 b"AC".to_vec(),
//                 b"AF".to_vec(),
//                 b"AN".to_vec(),
//                 b"DP".to_vec()
//             ]
//         );
//         assert_eq!(
//             config.format_list,
//             vec![b"AD".to_vec(), b"DP".to_vec(), b"GT".to_vec()]
//         );
//         assert!(config
//             .canonical_list
//             .as_ref()
//             .unwrap()
//             .contains(&b"ENST00000384061.1".to_vec()));
//         assert!(config
//             .canonical_list
//             .as_ref()
//             .unwrap()
//             .contains(&b"ENST00000533876.1_1".to_vec()));

//         let matches = VCF2CSV
//             .config_subcommand(App::new("test"))
//             .get_matches_from_safe(
//                 ["prog", "-m", "-i", "AC", "-f", "AD", "--decode-genotype"].iter(),
//             )
//             .unwrap();
//         let config = create_config(&vcf_reader.header(), &matches)?;
//         assert_eq!(
//             config,
//             VCF2CSVConfig {
//                 split_multi_allelic: true,
//                 decoded_genotype: true,
//                 canonical_list: None,
//                 info_list: vec![b"AC".to_vec(),],
//                 format_list: vec![b"AD".to_vec()],
//                 replace_sample_name: None,
//                 group_names: None,
//             }
//         );
//         Ok(())
//     }
// }
