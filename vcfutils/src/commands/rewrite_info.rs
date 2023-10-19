use crate::logic::rewrite_info::rewrite_info;
use crate::utils;
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
    // fn command_name(&self) -> &'static str {
    //     "rewrite-info"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("Rewrite INFO field")
    //         .arg(
    //             Arg::with_name("input")
    //                 .index(1)
    //                 .takes_value(true)
    //                 .help("Input VCF file"),
    //         )
    //         .arg(
    //             Arg::with_name("output")
    //                 .short("o")
    //                 .long("output")
    //                 .takes_value(true)
    //                 .help("Output file"),
    //         )
    //         .arg(
    //             Arg::with_name("info")
    //                 .short("i")
    //                 .long("info-list")
    //                 .takes_value(true)
    //                 .help("Include list of info tags")
    //                 .multiple(true),
    //         )
    //         .arg(
    //             Arg::with_name("exclude")
    //                 .short("e")
    //                 .long("exclude-info-list")
    //                 .takes_value(true)
    //                 .help("Exclude list of format tags")
    //                 .multiple(true),
    //         )
    // }

    pub fn run(&self) -> anyhow::Result<()> {
        let mut vcf_reader = utils::open_vcf_from_path(self.input.as_deref())?;
        let mut writer = autocompress::create_or_stdout(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?;
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
