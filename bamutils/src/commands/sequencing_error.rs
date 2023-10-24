use crate::logic::sequencing_error::{Regions, SequencingErrorProcessor};
use bio::io::fasta;
use clap::Args;
use rust_htslib::bam;
use std::{fs::File, str};

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
    #[arg(
        short = 'x',
        long,
        help = "maximum length of insertion/deletion to count",
        default_value = "3"
    )]
    max_indel_length: usize,
    #[arg(short, long, help = "Known SNP/INDELs VCF to exclude")]
    known_variants: Option<String>,
    #[arg(short, long, help = "Target regions BED")]
    region_bed: Option<String>,
}

impl SequencingError {
    pub fn run(&self) -> anyhow::Result<()> {
        let writer = autocompress::create_or_stdout(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?;
        let mut csv_writer = csv::WriterBuilder::new().flexible(true).from_writer(writer);

        let mut bam = bam::Reader::from_path(&self.bam)?;
        bam.set_reference(&self.reference)?;
        let mut reference_fasta = fasta::IndexedReader::from_file(&self.reference)?;
        let known_variants = self.known_variants.as_deref().map(|x| {
            rust_htslib::bcf::IndexedReader::from_path(x).expect("Failed to open known variants")
        });

        let regions = if let Some(bed) = self.region_bed.as_deref() {
            Regions::load_from_bed(File::open(bed)?)?
        } else {
            Regions::create_from_fasta(&mut reference_fasta)
        };

        let mut processor = SequencingErrorProcessor::new(
            self.max_indel_length,
            known_variants,
            reference_fasta,
            regions,
        );
        eprintln!("Counting...");
        processor.add_bam(&mut bam, self.min_mapq)?;
        let counter = processor.count();
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
        for (k, v) in counter.total_reference_base.iter() {
            csv_writer.write_record(&[
                "Total Reference Base".to_string(),
                str::from_utf8(&[*k])?.to_string(),
                "".to_string(),
                format!("{}", v),
            ])?;
        }
        for (k, v) in counter.total_reference_triplet.iter() {
            csv_writer.write_record(&[
                "Total Reference Triplet".to_string(),
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
        for (k, v) in counter.short_insertion.iter() {
            csv_writer.write_record(&[
                "Insertion".to_string(),
                "".to_string(),
                str::from_utf8(&k)?.to_string(),
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
        for (k, v) in counter.short_deletion.iter() {
            csv_writer.write_record(&[
                "Deletion".to_string(),
                "".to_string(),
                str::from_utf8(&k)?.to_string(),
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
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_run() -> anyhow::Result<()> {
        let cli = SequencingError {
            bam: "./testdata/demo1.cram".to_string(),
            reference: "./testdata/ref/MT.fa".to_string(),
            output: Some("../target/sequencing-error-count.csv".to_string()),
            min_mapq: 0,
            max_indel_length: 10,
            known_variants: None,
            region_bed: None,
        };
        cli.run()?;
        Ok(())
    }
}
