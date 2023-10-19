mod extend_bed;
mod merge_bed;

use anyhow::Result;
use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Commands {
    ExtendBed(extend_bed::ExtendBed),
    MergeBed(merge_bed::MergeBed),
}

impl Commands {
    pub fn run(&self) -> Result<()> {
        match self {
            Commands::ExtendBed(x) => x.run(),
            Commands::MergeBed(x) => x.run(),
        }
    }
}
