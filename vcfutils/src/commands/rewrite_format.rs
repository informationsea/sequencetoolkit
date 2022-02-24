use super::Command;
use crate::logic::rewrite_format::rewrite_format;
use crate::utils;
use clap::{App, Arg, ArgMatches};
use std::collections::HashSet;

pub struct RewriteFormat;

impl Command for RewriteFormat {
    fn command_name(&self) -> &'static str {
        "rewrite-format"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Rewrite FORMAT field")
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
                Arg::with_name("format")
                    .short("f")
                    .long("format-list")
                    .takes_value(true)
                    .help("Include list of format tags")
                    .multiple(true),
            )
            .arg(
                Arg::with_name("exclude")
                    .short("e")
                    .long("exclude-info-list")
                    .takes_value(true)
                    .help("Exclude list of format tags")
                    .multiple(true),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        let mut vcf_reader = utils::open_vcf_from_path(matches.value_of("input"))?;
        let mut writer = autocompress::create_or_stdout(
            matches.value_of("output"),
            autocompress::CompressionLevel::Default,
        )?;
        let blacklist = matches
            .values_of("exclude")
            .map(|x| x.map(|x| x.as_bytes().to_vec()).collect::<HashSet<_>>())
            .unwrap_or_default();
        let format_tags = matches
            .values_of("format")
            .map(|x| x.map(|x| x.as_bytes().to_vec()).collect::<Vec<_>>())
            .unwrap_or_else(|| vcf_reader.header().format_list().cloned().collect())
            .iter()
            .cloned()
            .filter(|x| !blacklist.contains(x))
            .collect::<HashSet<_>>();

        rewrite_format(&mut vcf_reader, &mut writer, &format_tags)?;

        Ok(())
    }
}
