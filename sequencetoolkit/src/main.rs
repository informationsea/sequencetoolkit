use clap::{Parser, Subcommand};
use std::env;

#[derive(Debug, Subcommand)]
pub enum Commands {
    BamUtils(bamutils::BamUtils),
    BedUtils(bedutils::BEDUtils),
    GeneAnnot(geneannot::GeneAnnot),
    VcfUtils(vcfutils::VCFUtils),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::BamUtils(x) => x.run(),
            Commands::BedUtils(x) => x.run(),
            Commands::GeneAnnot(x) => x.run(),
            Commands::VcfUtils(x) => x.run(),
        }
    }
}

#[derive(Debug, Parser)]
#[command(version, about = "sequence toolkit")]
pub struct Cli {
    #[arg(short = 'v', long = "verbose", action= clap::ArgAction::Count, help="verbose level")]
    verbose: u8,
    #[command(subcommand)]
    commands: Commands,
}

fn main() -> anyhow::Result<()> {
    let matches = Cli::parse();

    match matches.verbose {
        1 => env::set_var("RUST_LOG", "info"),
        2 => env::set_var("RUST_LOG", "debug"),
        3 => env::set_var("RUST_LOG", "trace"),
        _ => {
            if env::var("RUST_LOG").is_err() {
                env::set_var("RUST_LOG", "warn")
            }
        }
    }

    pretty_env_logger::init();

    matches.commands.run()?;

    Ok(())
}
