use super::Command;
use crate::logic::add_af;
use crate::utils;
use clap::{App, Arg, ArgMatches};
use std::collections::hash_map::RandomState;
use std::collections::HashMap;

pub struct AddAF;

impl Command for AddAF {
    fn command_name(&self) -> &'static str {
        "add-af"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Add or rewrite allele frequency, count, number and genotype count")
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
                    .help("Output CSV"),
            )
            .arg(
                Arg::with_name("category")
                    .short("c")
                    .long("category")
                    .help("category mapping file (csv or tsv)")
                    .requires("id")
                    .requires("value")
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("id")
                    .short("i")
                    .long("id")
                    .help("ID column name in category mapping file")
                    .takes_value(true)
                    .multiple(true),
            )
            .arg(
                Arg::with_name("value")
                    .short("v")
                    .long("value")
                    .help("value column name in category mapping file")
                    .takes_value(true)
                    .multiple(true),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        let mut vcf_reader = utils::open_vcf_from_path(matches.value_of("input"))?;
        let mut vcf_writer = autocompress::create_or_stdout(matches.value_of("output"))?;

        let category_to_sample = if let Some(x) = matches.value_of("category") {
            add_af::load_category_mapping::<_, RandomState>(
                &mut utils::auto_csv_reader_from_path(x, true)?,
                matches
                    .values_of("id")
                    .unwrap()
                    .map(|x| x.as_bytes().to_vec())
                    .collect(),
                matches
                    .values_of("value")
                    .unwrap()
                    .map(|x| x.as_bytes().to_vec())
                    .collect(),
            )?
        } else {
            HashMap::new()
        };

        add_af::add_af(&mut vcf_reader, &mut vcf_writer, &category_to_sample)?;

        Ok(())
    }
}
