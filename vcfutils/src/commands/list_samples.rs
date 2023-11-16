use crate::utils;
use autocompress::io::RayonWriter;
use clap::Args;
use std::io::prelude::*;
use std::str;

#[derive(Debug, Args)]
#[command(about = "List up sample names", version, author)]
pub struct ListSamples {
    #[arg(help = "Input VCF file")]
    input: Option<String>,
    #[arg(short, long, help = "Output Text")]
    output: Option<String>,
}

impl ListSamples {
    pub fn run(&self) -> anyhow::Result<()> {
        let vcf_reader = utils::open_vcf_from_path(self.input.as_deref())?;
        let mut output = RayonWriter::new(autocompress::autodetect_create_or_stdout_prefer_bgzip(
            self.output.as_ref(),
            autocompress::CompressionLevel::Default,
        )?);
        for one in vcf_reader.header().samples() {
            writeln!(output, "{}", str::from_utf8(one)?)?;
        }
        Ok(())
    }
}
