use crate::logic::sequencing_error::SequencingErrorCount;
use bio::io::fasta;
use clap::Args;
use rust_htslib::bam;
use std::str;

#[derive(Debug, Args)]
#[command(about = "Count sequencing error", version, author)]
pub struct SequencingError {
    #[arg(help = "Input BAM/CRAM file")]
    bam: String,
    #[arg(short = 'T', long, help = "Reference FASTA")]
    reference: String,
    #[arg(short, long, help = "Output CSV file")]
    output: Option<String>,
    #[arg(short, long, help = "Minimum mapping quality", default_value = "20")]
    min_mapq: u8,
}

impl SequencingError {
    pub fn run(&self) -> anyhow::Result<()> {
        run(
            &self.bam,
            &self.reference,
            self.output.as_deref(),
            self.min_mapq,
        )?;
        Ok(())
    }
}

fn run(
    bam_path: &str,
    reference_fasta_path: &str,
    output_path: Option<&str>,
    min_mapq: u8,
) -> anyhow::Result<()> {
    let writer =
        autocompress::create_or_stdout(output_path, autocompress::CompressionLevel::Default)?;
    let mut csv_writer = csv::WriterBuilder::new().flexible(true).from_writer(writer);

    let mut bam = bam::Reader::from_path(bam_path)?;
    bam.set_reference(reference_fasta_path)?;
    let mut reference_fasta = fasta::IndexedReader::from_file(&reference_fasta_path)?;
    let mut counter = SequencingErrorCount::new();
    eprintln!("Counting...");
    counter.add_bam(&mut bam, &mut reference_fasta, min_mapq)?;
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
            0,
        )?;
        Ok(())
    }
}
