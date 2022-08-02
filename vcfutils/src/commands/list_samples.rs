use super::Command;
use crate::utils;
use clap::{App, Arg, ArgMatches};
use std::str;

pub struct ListSamples;

impl Command for ListSamples {
    fn command_name(&self) -> &'static str {
        "list-samples"
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
                    .help("Output Text"),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        let vcf_reader = utils::open_vcf_from_path(matches.value_of("input"))?;
        let mut output = autocompress::create_or_stdout(
            matches.value_of("output"),
            autocompress::CompressionLevel::Default,
        )?;
        for one in vcf_reader.header().samples() {
            writeln!(output, "{}", str::from_utf8(one)?)?;
        }
        Ok(())
    }
}
