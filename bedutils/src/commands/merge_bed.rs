use super::super::logic::BedMerger;
use super::super::{BedReader, BedRegion};
use crate::Command;
use clap::{App, Arg, ArgMatches};
use std::io;
use std::str;

pub struct MergeBed;

impl Command for MergeBed {
    fn command_name(&self) -> &'static str {
        "merge"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Merge BED regions")
            .arg(
                Arg::with_name("input")
                    .index(1)
                    .takes_value(true)
                    .multiple(true)
                    .help("Input BED files"),
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
        let mut output = autocompress::create_or_stdout(
            matches.value_of("output"),
            autocompress::CompressionLevel::Default,
        )?;
        let mut merger = BedMerger::new();
        if let Some(inputs) = matches.values_of("input") {
            for one in inputs {
                let mut input = BedReader::new(io::BufReader::new(autocompress::open(one)?));
                let mut region = BedRegion::new();
                while input.next(&mut region)? {
                    merger.add(&region)
                }
            }
        } else {
            let stdin = io::stdin();
            let mut input = BedReader::new(io::BufReader::new(stdin.lock()));
            let mut region = BedRegion::new();
            while input.next(&mut region)? {
                merger.add(&region)
            }
        }
        merger.export_bed(&mut output)?;
        Ok(())
    }
}
