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

use crate::Command;

pub(crate) const COMMANDS: &[&dyn Command] = &[
    &vcf2csv::VCF2CSV,
    &add_af::AddAF,
    &list_samples::ListSamples,
    &rewrite_info::RewriteInfo,
    &rewrite_format::RewriteFormat,
    &remove_nonstandard_header::RemoveNonStandardHeader,
    &replace_contig::ReplaceContig,
    &generate_sql::GenerateSql,
    &replace_sample::ReplaceSampleName,
    &extract_canonical::ExtractCanonical,
];
