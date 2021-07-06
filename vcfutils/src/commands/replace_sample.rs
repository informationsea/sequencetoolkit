use super::Command;
use crate::error::VCFUtilsError;
use crate::logic::replace_sample::replace_sample;
use crate::utils::{self, Mapping};
use clap::{App, Arg, ArgMatches};
use rand::prelude::*;
use std::collections::HashMap;
use std::io::BufRead;
use vcf::{U8Vec, VCFReader};

pub struct ReplaceSampleName;

impl Command for ReplaceSampleName {
    fn command_name(&self) -> &'static str {
        "replace-sample-names"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Replace sample names")
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
                    .help("Output VCF file"),
            )
            .arg(
                Arg::with_name("mapping")
                    .short("m")
                    .long("sample-mapping")
                    .takes_value(true)
                    .help("sample name mapping file (csv/tsv)")
                    .long_help("sample name mapping file (csv/tsv). First column should be names to replace and second column should be new names. No header is required.")
                    .required_unless_one(&["sequential"]),
            )
            .arg(
                Arg::with_name("random")
                    .short("r")
                    .long("random")
                    .help("Shuffle sample order"),
            )
            .arg(
                Arg::with_name("sequential")
                    .short("s")
                    .long("rename-sequential")
                    .help("Rename sample as \"sample_1\", \"sample_2\", ..."),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        let mut reader = utils::open_vcf_from_path(matches.value_of("input"))?;
        let mut writer = autocompress::create_or_stdout(matches.value_of("output"))?;

        let (mapping, order): (HashMap<U8Vec, U8Vec>, Vec<U8Vec>) = create_sample_mapping(
            &reader,
            matches.value_of("mapping"),
            matches.is_present("sequential"),
            matches.is_present("random"),
        )?;
        replace_sample(&mut reader, &mut writer, &mapping, &order)?;
        Ok(())
    }
}

fn create_sample_mapping<R: BufRead>(
    reader: &VCFReader<R>,
    mapping_file: Option<&str>,
    sequential: bool,
    random: bool,
) -> Result<(HashMap<U8Vec, U8Vec>, Vec<U8Vec>), VCFUtilsError> {
    let mut rnd = thread_rng();
    Ok(if sequential {
        let mut original_names = reader.header().samples().to_vec();
        if random {
            original_names.shuffle(&mut rnd);
        }
        let new_names: Vec<_> = original_names
            .iter()
            .enumerate()
            .map(|(i, _)| format!("sample_{}", i + 1).into_bytes())
            .collect();
        (
            original_names
                .iter()
                .zip(new_names.iter())
                .map(|(x, y)| (x.clone(), y.clone()))
                .collect(),
            new_names,
        )
    } else {
        let Mapping {
            mapping,
            mut value_order,
            ..
        } = utils::load_mapping(utils::auto_csv_reader_from_path(
            mapping_file.unwrap(),
            false,
        )?)?;
        if random {
            value_order.shuffle(&mut rnd);
        }
        (mapping, value_order)
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_create_mapping() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/1kGP-subset.vcf");
        let vcf_reader = VCFReader::new(&vcf_data[..])?;
        let (mapping, order) = create_sample_mapping(&vcf_reader, None, true, false)?;
        let order_expected = vec![
            b"sample_1".to_vec(),
            b"sample_2".to_vec(),
            b"sample_3".to_vec(),
            b"sample_4".to_vec(),
            b"sample_5".to_vec(),
        ];
        assert_eq!(order, order_expected);

        let mapping_expected = [
            (b"SRP150637__HG00099", b"sample_1"),
            (b"SRP150637__HG00100", b"sample_2"),
            (b"SRP150637__HG00102", b"sample_3"),
            (b"SRP150637__HG00104", b"sample_4"),
            (b"SRP150637__HG00106", b"sample_5"),
        ]
        .iter()
        .map(|(x, y)| (x.to_vec(), y.to_vec()))
        .collect::<HashMap<_, _>>();
        assert_eq!(mapping, mapping_expected);

        let (mapping1, order) = create_sample_mapping(&vcf_reader, None, true, true)?;
        assert_eq!(order, order_expected);
        let (mapping2, order) = create_sample_mapping(&vcf_reader, None, true, true)?;
        assert_eq!(order, order_expected);
        let (mapping3, order) = create_sample_mapping(&vcf_reader, None, true, true)?;
        assert_eq!(order, order_expected);
        assert!(
            mapping_expected != mapping1
                || mapping_expected != mapping2
                || mapping_expected != mapping3
        );

        // with mapping file
        let (mapping, order) = create_sample_mapping(
            &vcf_reader,
            Some("./testfiles/1kGP-mapping.csv"),
            false,
            false,
        )?;

        let order_expected = vec![
            b"A".to_vec(),
            b"B".to_vec(),
            b"C".to_vec(),
            b"D".to_vec(),
            b"E".to_vec(),
        ];
        assert_eq!(order, order_expected);

        let mapping_expected = [
            (b"SRP150637__HG00099", b"A"),
            (b"SRP150637__HG00100", b"B"),
            (b"SRP150637__HG00102", b"C"),
            (b"SRP150637__HG00104", b"D"),
            (b"SRP150637__HG00106", b"E"),
        ]
        .iter()
        .map(|(x, y)| (x.to_vec(), y.to_vec()))
        .collect::<HashMap<_, _>>();
        assert_eq!(mapping, mapping_expected);

        Ok(())
    }
}
