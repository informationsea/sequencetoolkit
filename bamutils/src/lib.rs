mod commands;
pub mod logic;

use clap::Args;

#[derive(Debug, Args)]
#[command(about = "BAM Utilities", version, author)]
pub struct BamUtils {
    #[command(subcommand)]
    command: commands::Commands,
}

impl BamUtils {
    // fn command_name(&self) -> &'static str {
    //     "bamutils"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("BAM Utilities")
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
        // for one_command in commands::COMMANDS {
        //     if let Some(matches) = matches.subcommand_matches(one_command.command_name()) {
        //         return one_command.run(matches);
        //     }
        // }
        // unreachable!()
    }
}
