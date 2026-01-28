use clap::Args;
use rand::prelude::*;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::str;

#[derive(Debug, Args)]
#[command(about = "Random sample read", version, author)]
pub struct RandomSampling {
    #[arg(help = "Input BAM/CRAM file")]
    bam: String,
    #[arg(short = 'r', long, help = "Sampling rate")]
    sampling_rate: f64,
    #[arg(short = 'T', long, help = "Reference FASTA")]
    reference: Option<String>,
    #[arg(short, long, help = "Output BAM/CRAM/SAM file")]
    output: String,
    #[arg(short, long, help = "# of threads to write output file")]
    threads: Option<usize>,
}

impl RandomSampling {
    pub fn run(&self) -> anyhow::Result<()> {
        run(
            &self.bam,
            self.sampling_rate,
            &self.output,
            self.reference.as_deref(),
            self.threads,
        )?;
        Ok(())
    }
}

fn run(
    bam_path: &str,
    sampling_rate: f64,
    output_path: &str,
    reference: Option<&str>,
    threads: Option<usize>,
) -> anyhow::Result<()> {
    if sampling_rate <= 0. {
        return Err(anyhow::anyhow!("sampling rate must be larger than 0").into());
    }
    if sampling_rate > 1. {
        return Err(anyhow::anyhow!("sampling rate must be smaller than 1").into());
    }

    let mut rng = rand::rng();
    let mut reader = bam::Reader::from_path(bam_path)?;
    if let Some(reference) = reference {
        reader.set_reference(reference)?;
    }

    let output_format = if output_path.ends_with(".cram") {
        bam::Format::Cram
    } else if output_path.ends_with(".bam") {
        bam::Format::Bam
    } else if output_path.ends_with(".sam") {
        bam::Format::Sam
    } else {
        return Err(anyhow::anyhow!("Failed to detect output format: {}", output_path).into());
    };
    let mut output_header = bam::Header::from_template(reader.header());
    output_header.push_comment(b"random sampling with sequence toolkit");
    let mut writer = bam::Writer::from_path(output_path, &output_header, output_format)?;
    if let Some(threads) = threads {
        writer.set_threads(threads)?;
    }
    if let Some(reference) = reference {
        writer.set_reference(reference)?;
    }

    let mut decision_result: HashMap<Vec<u8>, bool> = HashMap::new();

    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        if let Some(d) = decision_result.get(record.qname()) {
            if *d {
                writer.write(&record)?;
            }
        } else {
            let r: f64 = rng.random();
            if r < sampling_rate {
                writer.write(&record)?;
                decision_result.insert(record.qname().to_vec(), true);
            } else {
                decision_result.insert(record.qname().to_vec(), false);
            }
        }
    }

    Ok(())
}
