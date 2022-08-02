use clap::{App, Arg, ArgMatches};
use rand::prelude::*;
use rust_htslib::bam::{self, Read};
use sequencetoolkit_common::Command;
use std::collections::HashMap;
use std::str;

pub struct RandomSampling;

impl Command for RandomSampling {
    fn command_name(&self) -> &'static str {
        return "random-sampling";
    }

    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Count sequencing error")
            .arg(
                Arg::with_name("bam")
                    .index(1)
                    .takes_value(true)
                    .required(true)
                    .help("Input BAM/CRAM file"),
            )
            .arg(
                Arg::with_name("sampling-rate")
                    .short("r")
                    .long("sampling-rate")
                    .required(true)
                    .takes_value(true)
                    .help("Sampling rate"),
            )
            .arg(
                Arg::with_name("reference")
                    .short("T")
                    .long("reference")
                    .takes_value(true)
                    .help("Reference FASTA"),
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("output")
                    .required(true)
                    .takes_value(true),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        run(
            matches.value_of("bam").unwrap(),
            matches
                .value_of("sampling-rate")
                .unwrap()
                .parse()
                .expect("Sampling rate must be float number"),
            matches.value_of("output").unwrap(),
            matches.value_of("reference"),
        )?;
        Ok(())
    }
}

fn run(
    bam_path: &str,
    sampling_rate: f64,
    output_path: &str,
    reference: Option<&str>,
) -> anyhow::Result<()> {
    if sampling_rate <= 0. {
        return Err(anyhow::anyhow!("sampling rate must be larger than 0").into());
    }
    if sampling_rate > 1. {
        return Err(anyhow::anyhow!("sampling rate must be smaller than 1").into());
    }

    let mut rng = rand::thread_rng();
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

    let mut decision_result: HashMap<Vec<u8>, bool> = HashMap::new();

    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        if let Some(d) = decision_result.get(record.qname()) {
            if *d {
                writer.write(&record)?;
            }
        } else {
            let r: f64 = rng.gen();
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
