mod fastq_triplet;

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Commands {
    FastqTriplet(fastq_triplet::FastqTriplet),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::FastqTriplet(x) => x.run(),
        }
    }
}
