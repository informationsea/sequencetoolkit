mod commands;
pub mod logic;

use clap::{crate_authors, crate_version, App, AppSettings, ArgMatches};
use std::str;

use sequencetoolkit_common::Command;

pub struct BamUtils;

impl Command for BamUtils {
    fn command_name(&self) -> &'static str {
        "bamutils"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("BAM Utilities")
            .version(crate_version!())
            .author(crate_authors!())
            .subcommands(commands::COMMANDS.iter().map(|x| {
                x.cli()
                    .setting(AppSettings::ColorAuto)
                    .setting(AppSettings::ColoredHelp)
            }))
            .setting(AppSettings::SubcommandRequiredElseHelp)
    }

    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        for one_command in commands::COMMANDS {
            if let Some(matches) = matches.subcommand_matches(one_command.command_name()) {
                return one_command.run(matches);
            }
        }
        unreachable!()
    }
}
