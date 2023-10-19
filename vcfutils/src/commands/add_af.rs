use crate::logic::add_af;
use crate::utils;
use anyhow::Context;
use clap::Args;
use std::collections::hash_map::RandomState;
use std::collections::HashMap;

#[derive(Debug, Args)]
#[command(
    about = "Add or rewrite allele frequency, count, number and genotype count",
    version,
    author
)]
pub struct AddAF {
    #[arg(help = "Input VCF file")]
    input: Option<String>,
    #[arg(short, long, help = "Output VCF file")]
    output: Option<String>,
    #[arg(
        short = 'c',
        long,
        help = "category mapping file (csv or tsv)",
        requires_all = &["id", "value"]
    )]
    category: Option<String>,
    #[arg(
        short = 'i',
        long,
        help = "ID column name in category mapping file",
        requires = "category"
    )]
    id: Option<Vec<String>>,
    #[arg(
        short = 'v',
        long,
        help = "value column name in category mapping file",
        requires = "category"
    )]
    value: Option<Vec<String>>,
    #[arg(
        short = 'p',
        long,
        help = "Allele frequency precision",
        default_value = "4"
    )]
    precision: usize,
}

impl AddAF {
    pub fn run(&self) -> anyhow::Result<()> {
        let mut vcf_reader = utils::open_vcf_from_path(self.input.as_deref())?;
        let mut vcf_writer = autocompress::create_or_stdout(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?;
        let af_precision: usize = self.precision;

        let category_to_sample = if let Some(x) = self.category.as_ref() {
            add_af::load_category_mapping::<_, RandomState>(
                &mut utils::auto_csv_reader_from_path(x, true)?,
                self.id
                    .as_ref()
                    .map(|x| x.iter().map(|x| x.as_bytes().to_vec()).collect())
                    .context("No ID column name")?,
                self.value
                    .as_ref()
                    .map(|x| x.iter().map(|x| x.as_bytes().to_vec()).collect())
                    .context("No value column name")?,
            )?
        } else {
            HashMap::new()
        };

        add_af::add_af(
            &mut vcf_reader,
            &mut vcf_writer,
            &category_to_sample,
            af_precision,
        )?;

        Ok(())
    }
}
