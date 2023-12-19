mod abspath;
mod autocat;
mod filterhash;

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[command(name = "abspath", alias = "ap")]
    AbsolutePath(abspath::AbsolutePath),
    #[command(name = "autocat", alias = "ac")]
    AutoCat(autocat::AutoCat),
    #[command(name = "filterhash", alias = "fh")]
    FilterHash(filterhash::FilterHash),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::AbsolutePath(x) => x.run(),
            Commands::AutoCat(x) => x.run(),
            Commands::FilterHash(x) => x.run(),
        }
    }
}
