use crate::logic::add_contig;
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
        let vcf_reader = autocompress::open(&self.input)?;
        let mut vcf_writer = autocompress::create_or_stdout(
            self.output.as_ref(),
            autocompress::CompressionLevel::Default,
        )?;

        let contigs = add_contig::scan_contig(&mut BufReader::new(vcf_reader))?;
        let vcf_reader = autocompress::open(&self.input)?;
        add_contig::add_contig(&mut BufReader::new(vcf_reader), &mut vcf_writer, &contigs)?;

        Ok(())
    }
}
