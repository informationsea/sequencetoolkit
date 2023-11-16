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
    pub fn run(&self) -> anyhow::Result<()> {
        self.commands.run()
    }
}
