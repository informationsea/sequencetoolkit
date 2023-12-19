use anyhow::Context;
use autocompress::autodetect_open;
use clap::Parser;

#[derive(Parser, Debug, Clone, PartialEq)]
#[command(
    about = "Automatically detect input file type, decompress and print it to stdout",
    version,
    author
)]
pub struct AutoCat {
    #[arg(help = "Path")]
    path: Vec<String>,
    #[arg(short, long, help = "Output path")]
    output: Option<String>,
}

impl AutoCat {
    pub fn run(&self) -> anyhow::Result<()> {
        let mut output_writer = autocompress::autodetect_create_or_stdout(
            self.output.as_ref(),
            autocompress::CompressionLevel::Default,
        )
        .with_context(|| {
            format!(
                "Failed to create {}",
                self.output.as_deref().unwrap_or("/dev/stdin")
            )
        })?;

        if self.path.is_empty() {
            let mut reader = autocompress::autodetect_reader(std::io::stdin().lock())?;
            std::io::copy(&mut reader, &mut output_writer)?;
        } else {
            for path in &self.path {
                if path == "-" {
                    let mut reader = autocompress::autodetect_reader(std::io::stdin().lock())?;
                    std::io::copy(&mut reader, &mut output_writer)?;
                } else {
                    let mut reader = autodetect_open(path)
                        .with_context(|| format!("Failed to open {}", path))?;
                    std::io::copy(&mut reader, &mut output_writer)?;
                }
            }
        }
        Ok(())
    }
}
