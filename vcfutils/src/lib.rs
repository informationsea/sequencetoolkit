use clap::{crate_authors, crate_version, App, AppSettings, ArgMatches};

pub mod commands;
pub mod error;
pub mod logic;
pub mod utils;

use sequencetoolkit_common::Command;
use sequencetoolkit_common::SequenceToolkitError;

pub struct VCFUtils;

impl Command for VCFUtils {
    fn command_name(&self) -> &'static str {
        "vcfutils"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("VCF Utilities")
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
