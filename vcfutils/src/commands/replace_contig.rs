use crate::logic::replace_contig::replace_contig;
use crate::utils;
use anyhow::Context;
use autocompress::io::{RayonReader, RayonWriter};
use clap::Args;
use std::collections::HashMap;
use vcf::U8Vec;

#[derive(Debug, Args)]
#[command(about = "Replace contig", version, author)]
pub struct ReplaceContig {
    #[arg(help = "Input VCF file")]
    input: Option<String>,
    #[arg(short, long, help = "Output VCF")]
    output: Option<String>,
    #[arg(short, long, help = "contig mapping file (csv/tsv)")]
    contig_mapping: Option<String>,
    #[arg(
        short = 'e',
        long,
        help = "Replace RefSeq contig name with chromosome number (for dbSNP 152/153)"
    )]
    replace_refseq: bool,
    #[arg(short, long, help = "add \"chr\" prefix to human chromosomes")]
    add_chr_prefix: bool,
    #[arg(short, long, help = "remove \"chr\" prefix to human chromosomes")]
    remove_chr_prefix: bool,
}

impl ReplaceContig {
    pub fn run(&self) -> anyhow::Result<()> {
        let mut reader =
            RayonReader::new(autocompress::autodetect_open_or_stdin(self.input.clone())?);
        let mut writer = RayonWriter::new(autocompress::autodetect_create_or_stdout_prefer_bgzip(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?);
        let mapping = if self.replace_refseq {
            refseq_replace()
        } else if self.add_chr_prefix {
            add_chr()
        } else if self.remove_chr_prefix {
            remove_chr()
        } else {
            utils::load_mapping(utils::auto_csv_reader_from_path(
                self.contig_mapping
                    .as_deref()
                    .context("No contig mapping file")?,
                false,
            )?)?
            .mapping
        };

        replace_contig(&mut reader, &mut writer, &mapping)?;
        Ok(())
    }
}

fn add_chr() -> HashMap<U8Vec, U8Vec> {
    [
        (b"1".to_vec(), b"chr1".to_vec()),
        (b"2".to_vec(), b"chr2".to_vec()),
        (b"3".to_vec(), b"chr3".to_vec()),
        (b"4".to_vec(), b"chr4".to_vec()),
        (b"5".to_vec(), b"chr5".to_vec()),
        (b"6".to_vec(), b"chr6".to_vec()),
        (b"7".to_vec(), b"chr7".to_vec()),
        (b"8".to_vec(), b"chr8".to_vec()),
        (b"9".to_vec(), b"chr9".to_vec()),
        (b"10".to_vec(), b"chr10".to_vec()),
        (b"11".to_vec(), b"chr11".to_vec()),
        (b"12".to_vec(), b"chr12".to_vec()),
        (b"13".to_vec(), b"chr13".to_vec()),
        (b"14".to_vec(), b"chr14".to_vec()),
        (b"15".to_vec(), b"chr15".to_vec()),
        (b"16".to_vec(), b"chr16".to_vec()),
        (b"17".to_vec(), b"chr17".to_vec()),
        (b"18".to_vec(), b"chr18".to_vec()),
        (b"19".to_vec(), b"chr19".to_vec()),
        (b"20".to_vec(), b"chr20".to_vec()),
        (b"21".to_vec(), b"chr21".to_vec()),
        (b"22".to_vec(), b"chr22".to_vec()),
        (b"X".to_vec(), b"chrX".to_vec()),
        (b"Y".to_vec(), b"chrY".to_vec()),
    ]
    .iter()
    .cloned()
    .collect()
}

fn remove_chr() -> HashMap<U8Vec, U8Vec> {
    [
        (b"chr1".to_vec(), b"1".to_vec()),
        (b"chr2".to_vec(), b"2".to_vec()),
        (b"chr3".to_vec(), b"3".to_vec()),
        (b"chr4".to_vec(), b"4".to_vec()),
        (b"chr5".to_vec(), b"5".to_vec()),
        (b"chr6".to_vec(), b"6".to_vec()),
        (b"chr7".to_vec(), b"7".to_vec()),
        (b"chr8".to_vec(), b"8".to_vec()),
        (b"chr9".to_vec(), b"9".to_vec()),
        (b"chr10".to_vec(), b"10".to_vec()),
        (b"chr11".to_vec(), b"11".to_vec()),
        (b"chr12".to_vec(), b"12".to_vec()),
        (b"chr13".to_vec(), b"13".to_vec()),
        (b"chr14".to_vec(), b"14".to_vec()),
        (b"chr15".to_vec(), b"15".to_vec()),
        (b"chr16".to_vec(), b"16".to_vec()),
        (b"chr17".to_vec(), b"17".to_vec()),
        (b"chr18".to_vec(), b"18".to_vec()),
        (b"chr19".to_vec(), b"19".to_vec()),
        (b"chr20".to_vec(), b"20".to_vec()),
        (b"chr21".to_vec(), b"21".to_vec()),
        (b"chr22".to_vec(), b"22".to_vec()),
        (b"chrX".to_vec(), b"X".to_vec()),
        (b"chrY".to_vec(), b"Y".to_vec()),
    ]
    .iter()
    .cloned()
    .collect()
}

fn refseq_replace() -> HashMap<U8Vec, U8Vec> {
    [
        (b"NC_000001.10".to_vec(), b"1".to_vec()),
        (b"NC_000002.11".to_vec(), b"2".to_vec()),
        (b"NC_000003.11".to_vec(), b"3".to_vec()),
        (b"NC_000004.11".to_vec(), b"4".to_vec()),
        (b"NC_000005.9".to_vec(), b"5".to_vec()),
        (b"NC_000006.11".to_vec(), b"6".to_vec()),
        (b"NC_000007.13".to_vec(), b"7".to_vec()),
        (b"NC_000008.10".to_vec(), b"8".to_vec()),
        (b"NC_000009.11".to_vec(), b"9".to_vec()),
        (b"NC_000010.10".to_vec(), b"10".to_vec()),
        (b"NC_000011.9".to_vec(), b"11".to_vec()),
        (b"NC_000012.11".to_vec(), b"12".to_vec()),
        (b"NC_000013.10".to_vec(), b"13".to_vec()),
        (b"NC_000014.8".to_vec(), b"14".to_vec()),
        (b"NC_000015.9".to_vec(), b"15".to_vec()),
        (b"NC_000016.9".to_vec(), b"16".to_vec()),
        (b"NC_000017.10".to_vec(), b"17".to_vec()),
        (b"NC_000018.9".to_vec(), b"18".to_vec()),
        (b"NC_000019.9".to_vec(), b"19".to_vec()),
        (b"NC_000020.10".to_vec(), b"20".to_vec()),
        (b"NC_000021.8".to_vec(), b"21".to_vec()),
        (b"NC_000022.10".to_vec(), b"22".to_vec()),
        (b"NC_000023.10".to_vec(), b"X".to_vec()),
        (b"NC_000024.9".to_vec(), b"Y".to_vec()),
        (b"NC_000001.11".to_vec(), b"1".to_vec()),
        (b"NC_000002.12".to_vec(), b"2".to_vec()),
        (b"NC_000003.12".to_vec(), b"3".to_vec()),
        (b"NC_000004.12".to_vec(), b"4".to_vec()),
        (b"NC_000005.10".to_vec(), b"5".to_vec()),
        (b"NC_000006.12".to_vec(), b"6".to_vec()),
        (b"NC_000007.14".to_vec(), b"7".to_vec()),
        (b"NC_000008.11".to_vec(), b"8".to_vec()),
        (b"NC_000009.12".to_vec(), b"9".to_vec()),
        (b"NC_000010.11".to_vec(), b"10".to_vec()),
        (b"NC_000011.10".to_vec(), b"11".to_vec()),
        (b"NC_000012.12".to_vec(), b"12".to_vec()),
        (b"NC_000013.11".to_vec(), b"13".to_vec()),
        (b"NC_000014.9".to_vec(), b"14".to_vec()),
        (b"NC_000015.10".to_vec(), b"15".to_vec()),
        (b"NC_000016.10".to_vec(), b"16".to_vec()),
        (b"NC_000017.11".to_vec(), b"17".to_vec()),
        (b"NC_000018.10".to_vec(), b"18".to_vec()),
        (b"NC_000019.10".to_vec(), b"19".to_vec()),
        (b"NC_000020.11".to_vec(), b"20".to_vec()),
        (b"NC_000021.9".to_vec(), b"21".to_vec()),
        (b"NC_000022.11".to_vec(), b"22".to_vec()),
        (b"NC_000023.11".to_vec(), b"X".to_vec()),
        (b"NC_000024.10".to_vec(), b"Y".to_vec()),
    ]
    .iter()
    .cloned()
    .collect()
}
