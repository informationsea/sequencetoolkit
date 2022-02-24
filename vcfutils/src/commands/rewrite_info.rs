use super::Command;
use crate::logic::rewrite_info::rewrite_info;
use crate::utils;
use clap::{App, Arg, ArgMatches};
use std::collections::HashSet;

pub struct RewriteInfo;

impl Command for RewriteInfo {
    fn command_name(&self) -> &'static str {
        "rewrite-info"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Rewrite INFO field")
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
                Arg::with_name("info")
                    .short("i")
                    .long("info-list")
                    .takes_value(true)
                    .help("Include list of info tags")
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
        let info_tags = matches
            .values_of("info")
            .map(|x| x.map(|x| x.as_bytes().to_vec()).collect::<Vec<_>>())
            .unwrap_or_else(|| vcf_reader.header().info_list().cloned().collect())
            .iter()
            .cloned()
            .filter(|x| !blacklist.contains(x))
            .collect::<Vec<_>>();

        rewrite_info(&mut vcf_reader, &mut writer, &info_tags)?;

        Ok(())
    }
}
