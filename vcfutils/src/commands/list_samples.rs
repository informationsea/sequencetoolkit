use crate::utils;
use clap::Args;
use std::str;

#[derive(Debug, Args)]
#[command(about = "List up sample names", version, author)]
pub struct ListSamples {
    #[arg(help = "Input VCF file")]
    input: Option<String>,
    #[arg(short, long, help = "Output Text")]
    output: Option<String>,
}

impl ListSamples {
    // fn command_name(&self) -> &'static str {
    //     "list-samples"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("List up sample names")
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
    //                 .help("Output Text"),
    //         )
    // }

    pub fn run(&self) -> anyhow::Result<()> {
        let vcf_reader = utils::open_vcf_from_path(self.input.as_deref())?;
        let mut output = autocompress::create_or_stdout(
            self.output.as_ref(),
            autocompress::CompressionLevel::Default,
        )?;
        for one in vcf_reader.header().samples() {
            writeln!(output, "{}", str::from_utf8(one)?)?;
        }
        Ok(())
    }
}
