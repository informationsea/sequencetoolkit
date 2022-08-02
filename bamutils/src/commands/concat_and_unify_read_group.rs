use clap::{App, Arg, ArgMatches};
use once_cell::sync::Lazy;
use regex::Regex;
use rust_htslib::bam::{self, Read};
use sequencetoolkit_common::Command;
use std::io::{prelude::*, BufReader};
use std::str;

static ID_TAG: Lazy<Regex> = Lazy::new(|| Regex::new(r"\bID:([^\s]+)\b").unwrap());

pub struct ConcatAndUnifyReadGroup;

impl Command for ConcatAndUnifyReadGroup {
    fn command_name(&self) -> &'static str {
        return "cat";
    }

    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Concatenate BAMs and unify read group")
            .arg(
                Arg::with_name("inputs")
                    .index(1)
                    .takes_value(true)
                    .multiple(true)
                    .required(true)
                    .help("Input BAM/CRAM file"),
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
            .arg(
                Arg::with_name("read-group")
                    .short("g")
                    .long("read-group")
                    .takes_value(true)
                    .help("Read group tag")
                    .default_value(
                        "@RG\tID:SAMPLE\tSM:SAMPLE\tLB:SAMPLE\tPL:mixed\tPU:SAMPLE\tCN:mixed",
                    ),
            )
            .arg(
                Arg::with_name("threads")
                    .short("t")
                    .long("threads")
                    .takes_value(true)
                    .help("# of threads to write output file"),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        run(
            matches
                .values_of("inputs")
                .unwrap()
                .map(|x| x.to_string())
                .collect(),
            matches.value_of("output").unwrap(),
            matches.value_of("reference"),
            matches.value_of("read-group").unwrap(),
            matches
                .value_of("threads")
                .map(|x| x.parse().expect("number of threads must be number")),
        )?;
        Ok(())
    }
}

fn run(
    bam_paths: Vec<String>,
    output_path: &str,
    reference: Option<&str>,
    read_group: &str,
    threads: Option<usize>,
) -> anyhow::Result<()> {
    let mut reader = bam::Reader::from_path(&bam_paths[0])?;
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

    let read_group_id = ID_TAG
        .captures(read_group)
        .expect("No ID is found in read group")
        .get(1)
        .unwrap()
        .as_str();

    let mut header = bam::header::Header::new();
    header.push_record(&bam::header::HeaderRecord::new(b"HD\tVN:1.6"));

    let mut template_header = BufReader::new(reader.header().as_bytes());
    let mut header_line = Vec::new();
    while template_header.read_until(b'\n', &mut header_line)? > 0 {
        let line = header_line.strip_suffix(b"\n").unwrap_or(&header_line);
        if line.starts_with(b"@SQ") {
            header.push_record(&bam::header::HeaderRecord::new(&line[1..]));
        }
        header_line.clear();
    }
    header.push_record(&bam::header::HeaderRecord::new(&read_group.as_bytes()[1..]));
    header.push_record(&bam::header::HeaderRecord::new(b"PG\tID:sequencetoolkit\tPN:sequencetoolkit\tDS:concatenate bam files and unify read group."));

    let mut writer = bam::Writer::from_path(output_path, &header, output_format)?;
    if let Some(threads) = threads {
        writer.set_threads(threads)?;
    }

    let mut record = bam::Record::new();
    while let Some(r) = reader.read(&mut record) {
        r?;

        let mut found_rg = false;
        for x in record.aux_iter() {
            let x = x?;
            if x.0 == b"RG" {
                found_rg = true;
            }
        }
        if found_rg {
            record.remove_aux(b"RG")?;
        }

        record.push_aux(b"RG", bam::record::Aux::String(read_group_id))?;
        writer.write(&record)?;
    }
    std::mem::drop(reader);

    for one in &bam_paths[1..] {
        let mut reader = bam::Reader::from_path(one)?;
        while let Some(r) = reader.read(&mut record) {
            r?;
            writer.write(&record)?;
        }
    }

    Ok(())
}