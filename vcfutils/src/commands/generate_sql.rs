use super::vcf2csv::create_config;
use super::Command;
use crate::logic::generate_sql::generate_sql;
use crate::utils;
use clap::{App, Arg, ArgMatches};
use std::str;

pub struct GenerateSql;

impl Command for GenerateSql {
    fn command_name(&self) -> &'static str {
        "generate-sql"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("List up sample names")
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
                    .help("Output SQL"),
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
            .arg(
                Arg::with_name("table-name")
                    .short("n")
                    .long("table-name")
                    .help("Table name")
                    .takes_value(true)
                    .required(true),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        let vcf_reader = utils::open_vcf_from_path(matches.value_of("input"))?;
        let mut output = autocompress::create_or_stdout(
            matches.value_of("output"),
            autocompress::CompressionLevel::Default,
        )?;
        let mut config = create_config(&vcf_reader.header(), matches)?;
        config.split_multi_allelic = true;
        config.decoded_genotype = false;
        let generated_sql = generate_sql(
            &vcf_reader.header(),
            &config,
            matches.value_of("table-name").unwrap(),
        )?;
        output.write_all(generated_sql.as_bytes())?;
        Ok(())
    }
}
