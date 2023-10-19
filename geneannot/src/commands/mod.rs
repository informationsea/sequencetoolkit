mod create_db;
mod genome_position;
mod transcript_position;

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Commands {
    CreateDb(create_db::CreateDb),
    GenomePosition(genome_position::GenomePosition),
    TranscriptPosition(transcript_position::TranscriptPosition),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::CreateDb(x) => x.run(),
            Commands::GenomePosition(x) => x.run(),
            Commands::TranscriptPosition(x) => x.run(),
        }
    }
}
