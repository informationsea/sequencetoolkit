use clap::{crate_authors, crate_version, App, AppSettings, Arg};

use sequencetoolkit_common::Command;
use std::env;

pub(crate) const COMMANDS: &[&dyn Command] = &[
    &vcfutils::VCFUtils,
    &geneannot::GeneAnnot,
    &bedutils::BEDUtils,
];

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
