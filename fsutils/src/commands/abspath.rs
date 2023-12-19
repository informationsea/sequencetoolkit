use anyhow::Context;
use clap::Parser;

#[derive(Parser, Debug, Clone, PartialEq)]
#[command(about = "Print absolute path", version, author)]
pub struct AbsolutePath {
    #[arg(help = "Path")]
    path: Vec<String>,
}

impl AbsolutePath {
    pub fn run(&self) -> anyhow::Result<()> {
        for path in &self.path {
            println!(
                "{}",
                std::fs::canonicalize(path)
                    .context("Failed to get absolute path")?
                    .display()
            );
        }
        Ok(())
    }
}
