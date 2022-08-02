use crate::logic::sequencing_error::SequencingErrorCount;
use bio::io::fasta;
use clap::{App, Arg, ArgMatches};
use rust_htslib::bam;
use sequencetoolkit_common::Command;
use std::str;

pub struct SequencingError;

impl Command for SequencingError {
    fn command_name(&self) -> &'static str {
        return "sequencing-error";
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
                Arg::with_name("reference")
                    .short("r")
                    .long("reference")
                    .alias("T")
                    .required(true)
                    .takes_value(true)
                    .help("Reference FASTA"),
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("output")
                    .takes_value(true),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        run(
            matches.value_of("bam").unwrap(),
            matches.value_of("reference").unwrap(),
            matches.value_of("output"),
        )?;
        Ok(())
    }
}

fn run(
    bam_path: &str,
    reference_fasta_path: &str,
    output_path: Option<&str>,
) -> anyhow::Result<()> {
    let writer =
        autocompress::create_or_stdout(output_path, autocompress::CompressionLevel::Default)?;
    let mut csv_writer = csv::WriterBuilder::new().flexible(true).from_writer(writer);

    let mut bam = bam::Reader::from_path(bam_path)?;
    bam.set_reference(reference_fasta_path)?;
    let mut reference_fasta = fasta::IndexedReader::from_file(&reference_fasta_path)?;
    let mut counter = SequencingErrorCount::new();
    eprintln!("Counting...");
    counter.add_bam(&mut bam, &mut reference_fasta)?;
    eprintln!("Writing result...");
    csv_writer.write_record(&["Category", "Reference", "Sequenced", "Count"])?;
    csv_writer.write_record(&[
        "Total Sequenced Bases".to_string(),
        "".to_string(),
        "".to_string(),
        format!("{}", counter.total_sequenced_len),
    ])?;
    csv_writer.write_record(&[
        "Total Reference Bases".to_string(),
        "".to_string(),
        "".to_string(),
        format!("{}", counter.total_reference_len),
    ])?;
    for (k, v) in counter.total_match_base.iter() {
        csv_writer.write_record(&[
            "Total Match Reference Base".to_string(),
            str::from_utf8(&[*k])?.to_string(),
            "".to_string(),
            format!("{}", v),
        ])?;
    }
    for (k, v) in counter.total_match_triplet.iter() {
        csv_writer.write_record(&[
            "Total Match Reference Triplet".to_string(),
            str::from_utf8(k)?.to_string(),
            "".to_string(),
            format!("{}", v),
        ])?;
    }
    for (k, v) in counter.mismatch.iter() {
        csv_writer.write_record(&[
            "Mismatch".to_string(),
            str::from_utf8(&[k.reference])?.to_string(),
            str::from_utf8(&[k.sequenced])?.to_string(),
            format!("{}", v),
        ])?;
    }
    for (k, v) in counter.mismatch_triplet.iter() {
        csv_writer.write_record(&[
            "Mismatch Triplet".to_string(),
            str::from_utf8(&k.reference)?.to_string(),
            str::from_utf8(&k.sequenced)?.to_string(),
            format!("{}", v),
        ])?;
    }
    for (k, v) in counter.insertion_length.iter() {
        csv_writer.write_record(&[
            "Insertion Length".to_string(),
            "".to_string(),
            format!("{}", k),
            format!("{}", v),
        ])?;
    }
    for (k, v) in counter.deletion_length.iter() {
        csv_writer.write_record(&[
            "Deletion Length".to_string(),
            "".to_string(),
            format!("{}", k),
            format!("{}", v),
        ])?;
    }
    for (k, v) in counter.softclip_length.iter() {
        csv_writer.write_record(&[
            "Softclip Length".to_string(),
            "".to_string(),
            format!("{}", k),
            format!("{}", v),
        ])?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_run() -> anyhow::Result<()> {
        run(
            "./testdata/demo1.cram",
            "./testdata/ref/MT.fa",
            Some("../target/sequencing-error-count.csv"),
        )?;
        Ok(())
    }
}
