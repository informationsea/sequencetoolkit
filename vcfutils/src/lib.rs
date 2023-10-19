use clap::Args;

pub mod commands;
pub mod error;
pub mod logic;
pub mod utils;

#[derive(Debug, Args)]
#[command(about = "VCF Utilities", version, author)]
pub struct VCFUtils {
    #[command(subcommand)]
    commands: commands::Commands,
}

impl VCFUtils {
    // fn command_name(&self) -> &'static str {
    //     "vcfutils"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("VCF Utilities")
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
        self.commands.run()
    }
}
