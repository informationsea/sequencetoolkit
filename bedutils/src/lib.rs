mod commands;
pub mod logic;
mod parser;
pub mod utils;

use clap::Args;
pub use parser::{BedReader, BedRegion, BedWriter};

#[derive(Debug, Args)]
#[command(about = "BED Utilities", version, author)]
pub struct BEDUtils {
    #[command(subcommand)]
    command: commands::Commands,
}

impl BEDUtils {
    // fn command_name(&self) -> &'static str {
    //     "bedutils"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("BED Utilities")
    //         .version(crate_version!())
    //         .author(crate_authors!())
    //         .subcommands(commands::COMMANDS.iter().map(|x| {
    //             x.cli()
    //                 .setting(AppSettings::ColorAuto)
    //                 .setting(AppSettings::ColoredHelp)
    //         }))
    //         .setting(AppSettings::SubcommandRequiredElseHelp)
    // }

    pub fn run(&self) -> anyhow::Result<()> {
        self.command.run()
    }
}
