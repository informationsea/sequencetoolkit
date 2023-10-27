mod fastq_triplet;
mod index_count;

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Commands {
    FastqTriplet(fastq_triplet::FastqTriplet),
    IndexCount(index_count::IndexCount),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::FastqTriplet(x) => x.run(),
            Commands::IndexCount(x) => x.run(),
        }
    }
}
