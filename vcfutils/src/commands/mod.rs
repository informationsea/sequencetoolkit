mod add_af;
mod add_contig;
mod extract_canonical;
mod generate_sql;
mod list_samples;
mod remove_nonstandard_header;
mod replace_contig;
mod replace_sample;
mod rewrite_format;
mod rewrite_info;
mod vcf2csv;

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Commands {
    VCF2CSV(vcf2csv::VCF2CSV),
    AddAF(add_af::AddAF),
    AddContig(add_contig::AddContig),
    ListSamples(list_samples::ListSamples),
    RewriteInfo(rewrite_info::RewriteInfo),
    RewriteFormat(rewrite_format::RewriteFormat),
    RemoveNonStandardHeader(remove_nonstandard_header::RemoveNonStandardHeader),
    ReplaceContig(replace_contig::ReplaceContig),
    GenerateSql(generate_sql::GenerateSql),
    ReplaceSampleName(replace_sample::ReplaceSampleName),
    ExtractCanonical(extract_canonical::ExtractCanonical),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::VCF2CSV(x) => x.run(),
            Commands::AddAF(x) => x.run(),
            Commands::AddContig(x) => x.run(),
            Commands::ListSamples(x) => x.run(),
            Commands::RewriteInfo(x) => x.run(),
            Commands::RewriteFormat(x) => x.run(),
            Commands::RemoveNonStandardHeader(x) => x.run(),
            Commands::ReplaceContig(x) => x.run(),
            Commands::GenerateSql(x) => x.run(),
            Commands::ReplaceSampleName(x) => x.run(),
            Commands::ExtractCanonical(x) => x.run(),
        }
    }
}
