use clap::{crate_authors, crate_version, App, AppSettings, Arg, ArgMatches, SubCommand};

pub mod bedutils;
pub mod error;
pub mod geneannot;
pub mod vcfutils;

use error::{SequenceToolkitError, SequenceToolkitErrorKind};
use std::env;

pub(crate) const COMMANDS: &[&dyn Command] = &[
    &vcfutils::VCFUtils,
    &geneannot::GeneAnnot,
    &bedutils::BEDUtils,
];

pub trait Command {
    fn cli(&self) -> App<'static, 'static> {
        self.config_subcommand(SubCommand::with_name(self.command_name()))
            .version(crate_version!())
            .author(crate_authors!())
    }
    fn command_name(&self) -> &'static str;
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static>;
    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), SequenceToolkitError>;
}

fn main() {
    let matches = App::new("sequence toolkit")
        .version(crate_version!())
        .author(crate_authors!())
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .multiple(true),
        )
        .subcommands(COMMANDS.iter().map(|x| {
            x.cli()
                .setting(AppSettings::ColorAuto)
                .setting(AppSettings::ColoredHelp)
        }))
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::ColoredHelp)
        .get_matches();

    match matches.occurrences_of("verbose") {
        1 => env::set_var("RUST_LOG", "info"),
        2 => env::set_var("RUST_LOG", "debug"),
        3 => env::set_var("RUST_LOG", "trace"),
        _ => {
            if env::var("RUST_LOG").is_err() {
                env::set_var("RUST_LOG", "warn")
            }
        }
    }

    pretty_env_logger::init();

    for one_command in COMMANDS {
        if let Some(matches) = matches.subcommand_matches(one_command.command_name()) {
            one_command.run(matches).expect("Operation Error");
            return;
        }
    }
    unreachable!()
}
