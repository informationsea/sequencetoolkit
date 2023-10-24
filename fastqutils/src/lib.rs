pub mod commands;

use clap::Args;

#[derive(Debug, Args)]
#[command(version, about = "FASTQ Utilities", author)]
pub struct FastqUtils {
    #[command(subcommand)]
    commands: commands::Commands,
}

impl FastqUtils {
    pub fn run(&self) -> anyhow::Result<()> {
        self.commands.run()
    }
}
