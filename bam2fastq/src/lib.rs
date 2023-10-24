pub mod commands;
pub mod utils;

use clap::Args;

#[derive(Debug, Clone, Args)]
#[command(version, about = "Convert BAM/CRAM to FASTQ", author)]
pub struct Bam2fastq {
    #[command(subcommand)]
    command: commands::Commands,
    #[arg(short = 'v', long = "verbose", action= clap::ArgAction::Count, help="verbose level")]
    verbose: u8,
}

impl Bam2fastq {
    pub fn run(&self) -> anyhow::Result<()> {
        self.command.run()
    }
}
