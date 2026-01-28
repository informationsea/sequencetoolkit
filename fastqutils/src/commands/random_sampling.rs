use anyhow::Context;
use autocompress::io::{RayonReader, RayonWriter};
use clap::Parser;
use rand::Rng;
use std::io::prelude::*;

#[derive(Parser, Debug, Clone, PartialEq)]
#[command(about = "Random sampling of FASTQ files", version, author)]
pub struct RandomSampling {
    #[arg(short = '1', long, alias = "i1", help = "Input FASTQ 1 (with gzip)")]
    input1: String,
    #[arg(short = '2', long, help = "Input FASTQ 2 (with gzip)")]
    input2: String,
    #[arg(short = 'a', long, help = "Output FASTQ1 (with gzip)")]
    output1: String,
    #[arg(short = 'b', long, help = "Output FASTQ2 (with gzip)")]
    output2: String,
    #[arg(short, long, help = "Sampling ratio", default_value = "0.1")]
    ratio: f64,
    #[arg(short, long, help = "# of threads", default_value = "1")]
    threads: usize,
}

impl RandomSampling {
    pub fn run(&self) -> anyhow::Result<()> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .context("Failed to build thread pool")?;
        let mut rng = rand::rng();

        let mut reader1 = RayonReader::new(
            autocompress::autodetect_open(&self.input1)
                .with_context(|| format!("Failed to open {}", self.input1))?,
        );
        let mut reader2 = RayonReader::new(
            autocompress::autodetect_open(&self.input2)
                .with_context(|| format!("Failed to open {}", self.input2))?,
        );
        let mut writer1 = RayonWriter::new(
            autocompress::autodetect_create(
                &self.output1,
                autocompress::CompressionLevel::default(),
            )
            .with_context(|| format!("Failed to create {}", self.output1))?,
        );
        let mut writer2 = RayonWriter::new(
            autocompress::autodetect_create(
                &self.output2,
                autocompress::CompressionLevel::default(),
            )
            .with_context(|| format!("Failed to create {}", self.output2))?,
        );

        let mut line = String::new();
        loop {
            line.clear();
            if reader1.read_line(&mut line)? == 0 {
                if reader2.read_line(&mut line)? != 0 {
                    anyhow::bail!("FASTQ files have different number of lines");
                }
                break;
            }
            let rnd_value: f64 = rng.random();
            if rnd_value < self.ratio {
                writer1.write_all(line.as_bytes())?;
                for _ in 0..3 {
                    line.clear();
                    if reader1.read_line(&mut line)? == 0 {
                        anyhow::bail!("Unexpected end of FASTQ1")
                    }
                    writer1.write_all(line.as_bytes())?;
                }

                for _ in 0..4 {
                    line.clear();
                    if reader2.read_line(&mut line)? == 0 {
                        anyhow::bail!("Unexpected end of FASTQ2")
                    }
                    writer2.write_all(line.as_bytes())?;
                }
            } else {
                for _ in 0..3 {
                    line.clear();
                    if reader1.read_line(&mut line)? == 0 {
                        anyhow::bail!("Unexpected end of FASTQ1")
                    }
                }

                for _ in 0..4 {
                    line.clear();
                    if reader2.read_line(&mut line)? == 0 {
                        anyhow::bail!("Unexpected end of FASTQ2")
                    }
                }
            }
        }

        Ok(())
    }
}
