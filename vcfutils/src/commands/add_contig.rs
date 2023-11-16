use crate::logic::add_contig;
use autocompress::io::{RayonReader, RayonWriter};
use clap::Args;
use std::io::BufReader;

#[derive(Debug, Args)]
#[command(about = "Scan VCF file and add contig header", version, author)]
pub struct AddContig {
    #[arg(help = "Input VCF file")]
    input: String,
    #[arg(short, long, help = "Output VCF file")]
    output: Option<String>,
}

impl AddContig {
    pub fn run(&self) -> anyhow::Result<()> {
        let vcf_reader = RayonReader::new(autocompress::autodetect_open(&self.input)?);
        let mut vcf_writer =
            RayonWriter::new(autocompress::autodetect_create_or_stdout_prefer_bgzip(
                self.output.as_ref(),
                autocompress::CompressionLevel::Default,
            )?);

        let contigs = add_contig::scan_contig(&mut BufReader::new(vcf_reader))?;
        let vcf_reader = autocompress::autodetect_open(&self.input)?;
        add_contig::add_contig(&mut BufReader::new(vcf_reader), &mut vcf_writer, &contigs)?;

        Ok(())
    }
}
