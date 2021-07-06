pub mod annotator;
mod commands;
pub mod utils;

pub use annotator::error::{GeneAnnotError, GeneAnnotErrorKind};
use clap::{crate_authors, crate_version, App, AppSettings, ArgMatches};
use std::env;

use sequencetoolkit_common::{Command, SequenceToolkitError};

pub struct GeneAnnot;

impl Command for GeneAnnot {
    fn command_name(&self) -> &'static str {
        "geneannot"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Gene Annotator")
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
        eprintln!("subcommand is required");
        eprintln!("Run geneannot help to get help");
        Ok(())
    }
}
