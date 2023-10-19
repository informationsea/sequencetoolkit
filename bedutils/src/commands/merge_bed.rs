use super::super::logic::BedMerger;
use super::super::{BedReader, BedRegion};
use clap::Args;
use std::io;

#[derive(Debug, Args)]
pub struct MergeBed {
    #[arg(help = "Input BED files")]
    input: Option<Vec<String>>,
    #[arg(short, long, help = "Output BED file")]
    output: Option<String>,
}

impl MergeBed {
    // fn command_name(&self) -> &'static str {
    //     "merge"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("Merge BED regions")
    //         .arg(
    //             Arg::with_name("input")
    //                 .index(1)
    //                 .takes_value(true)
    //                 .multiple(true)
    //                 .help("Input BED files"),
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
        let mut output = autocompress::create_or_stdout(
            self.output.as_ref(),
            autocompress::CompressionLevel::Default,
        )?;
        let mut merger = BedMerger::new();
        if let Some(inputs) = self.input.as_ref() {
            for one in inputs {
                let mut input = BedReader::new(io::BufReader::new(autocompress::open(one)?));
                let mut region = BedRegion::new();
                while input.next(&mut region)? {
                    merger.add(&region)
                }
            }
        } else {
            let stdin = io::stdin();
            let mut input = BedReader::new(io::BufReader::new(stdin.lock()));
            let mut region = BedRegion::new();
            while input.next(&mut region)? {
                merger.add(&region)
            }
        }
        merger.export_bed(&mut output)?;
        Ok(())
    }
}
