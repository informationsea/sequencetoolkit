use crate::logic::rewrite_format::rewrite_format;
use crate::utils;
use clap::Args;
use std::collections::HashSet;

#[derive(Args, Debug)]
#[command(about = "Rewrite FORMAT field", version, author)]
pub struct RewriteFormat {
    #[arg(help = "Input VCF file")]
    input: Option<String>,
    #[arg(short, long, help = "Output VCF")]
    output: Option<String>,
    #[arg(short, long, help = "Include list of format tags")]
    format: Option<Vec<String>>,
    #[arg(short, long, help = "Exclude list of format tags")]
    exclude: Option<Vec<String>>,
}

impl RewriteFormat {
    // fn command_name(&self) -> &'static str {
    //     "rewrite-format"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("Rewrite FORMAT field")
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
    //             Arg::with_name("format")
    //                 .short("f")
    //                 .long("format-list")
    //                 .takes_value(true)
    //                 .help("Include list of format tags")
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
        let format_tags = self
            .format
            .as_ref()
            .map(|x| x.iter().map(|x| x.as_bytes().to_vec()).collect::<Vec<_>>())
            .unwrap_or_else(|| vcf_reader.header().format_list().cloned().collect())
            .iter()
            .cloned()
            .filter(|x| !blacklist.contains(x))
            .collect::<HashSet<_>>();

        rewrite_format(&mut vcf_reader, &mut writer, &format_tags)?;

        Ok(())
    }
}
