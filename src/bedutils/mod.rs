mod commands;
pub mod logic;
mod parser;
pub mod utils;

use crate::SequenceToolkitError;
use clap::{crate_authors, crate_version, App, AppSettings, ArgMatches};
use std::str;

pub use parser::{BedReader, BedRegion, BedWriter};

use crate::Command;

pub struct BEDUtils;

impl Command for BEDUtils {
    fn command_name(&self) -> &'static str {
        "bedutils"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("BED Utilities")
            .version(crate_version!())
            .author(crate_authors!())
            .subcommands(commands::COMMANDS.iter().map(|x| {
                x.cli()
                    .setting(AppSettings::ColorAuto)
                    .setting(AppSettings::ColoredHelp)
            }))
            .setting(AppSettings::SubcommandRequiredElseHelp)
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), SequenceToolkitError> {
        for one_command in commands::COMMANDS {
            if let Some(matches) = matches.subcommand_matches(one_command.command_name()) {
                return one_command.run(matches);
            }
        }
        unreachable!()
    }
}
