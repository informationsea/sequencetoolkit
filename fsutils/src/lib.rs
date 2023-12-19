pub mod commands;

use clap::Args;

#[derive(Debug, Args)]
#[command(version, about = "File System Utilities", author)]
pub struct FileSystemUtilities {
    #[command(subcommand)]
    commands: commands::Commands,
}

impl FileSystemUtilities {
    pub fn run(&self) -> anyhow::Result<()> {
        self.commands.run()
    }
}
