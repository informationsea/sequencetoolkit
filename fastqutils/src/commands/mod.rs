mod fastq_triplet;
mod index_count;
mod random_sampling;

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Commands {
    FastqTriplet(fastq_triplet::FastqTriplet),
    IndexCount(index_count::IndexCount),
    RandomSampling(random_sampling::RandomSampling),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::FastqTriplet(x) => x.run(),
            Commands::IndexCount(x) => x.run(),
            Commands::RandomSampling(x) => x.run(),
        }
    }
}
