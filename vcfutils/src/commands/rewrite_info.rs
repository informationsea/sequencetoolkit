use crate::logic::rewrite_info::rewrite_info;
use crate::utils;
use autocompress::io::RayonWriter;
use clap::Args;
use std::collections::HashSet;

#[derive(Args, Debug)]
#[command(about = "Replace INFO field", version, author)]
pub struct RewriteInfo {
    #[arg(help = "Input VCF file")]
    input: Option<String>,
    #[arg(short, long, help = "Output VCF")]
    output: Option<String>,
    #[arg(short, long, help = "Include list of info tags")]
    info: Option<Vec<String>>,
    #[arg(short, long, help = "Exclude list of format tags")]
    exclude: Option<Vec<String>>,
}

impl RewriteInfo {
    pub fn run(&self) -> anyhow::Result<()> {
        let mut vcf_reader = utils::open_vcf_from_path(self.input.as_deref())?;
        let mut writer = RayonWriter::new(autocompress::autodetect_create_or_stdout_prefer_bgzip(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?);
        let blacklist = self
            .exclude
            .as_ref()
            .map(|x| {
                x.iter()
                    .map(|x| x.as_bytes().to_vec())
                    .collect::<HashSet<_>>()
            })
            .unwrap_or_default();
        let info_tags = self
            .info
            .as_ref()
            .map(|x| x.iter().map(|x| x.as_bytes().to_vec()).collect::<Vec<_>>())
            .unwrap_or_else(|| vcf_reader.header().info_list().cloned().collect())
            .iter()
            .cloned()
            .filter(|x| !blacklist.contains(x))
            .collect::<Vec<_>>();

        rewrite_info(&mut vcf_reader, &mut writer, &info_tags)?;

        Ok(())
    }
}
