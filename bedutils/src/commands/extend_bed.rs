use super::super::{BedReader, BedRegion, BedWriter};
use clap::Args;
use std::io;
use std::u64;

#[derive(Debug, Args)]
pub struct ExtendBed {
    #[arg(help = "Input BED files")]
    input: Option<String>,
    #[arg(short, long, help = "Output BED")]
    output: Option<String>,
    #[arg(short = 'e', long, help = "Expand length", conflicts_with_all= ["expand_start", "expand_end"])]
    expand: Option<u64>,
    #[arg(
        short = 's',
        long,
        help = "Expand start length",
        requires = "expand_end"
    )]
    expand_start: Option<u64>,
    #[arg(
        short = 'n',
        long,
        help = "Expand end length",
        requires = "expand_start"
    )]
    expand_end: Option<u64>,
}

impl ExtendBed {
    // fn command_name(&self) -> &'static str {
    //     "expand"
    // }

    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("Extend BED regions")
    //         .arg(
    //             Arg::with_name("input")
    //                 .index(1)
    //                 .takes_value(true)
    //                 .help("Input BED files"),
    //         )
    //         .arg(
    //             Arg::with_name("output")
    //                 .short("o")
    //                 .long("output")
    //                 .takes_value(true)
    //                 .help("Output Text"),
    //         )
    //         .arg(
    //             Arg::with_name("expand")
    //                 .short("e")
    //                 .long("expand")
    //                 .takes_value(true)
    //                 .help("Expand length")
    //                 .required_unless("expand-start"),
    //         )
    //         .arg(
    //             Arg::with_name("expand-start")
    //                 .short("s")
    //                 .long("expand-start")
    //                 .takes_value(true)
    //                 .help("Expand start position")
    //                 .conflicts_with("expand")
    //                 .requires("expand-end"),
    //         )
    //         .arg(
    //             Arg::with_name("expand-end")
    //                 .short("n")
    //                 .long("expand-end")
    //                 .takes_value(true)
    //                 .help("Expand end position")
    //                 .conflicts_with("expand")
    //                 .requires("expand-start"),
    //         )
    // }

    pub fn run(&self) -> anyhow::Result<()> {
        let expand_start: u64 = self.expand_start.unwrap_or_else(|| self.expand.unwrap());
        let expand_end: u64 = self.expand_end.unwrap_or_else(|| self.expand.unwrap());

        let mut output = BedWriter::new(io::BufWriter::new(
            autocompress::autodetect_create_or_stdout(
                self.output.as_deref(),
                autocompress::CompressionLevel::Default,
            )?,
        ));
        let mut input = BedReader::new(io::BufReader::new(autocompress::autodetect_open_or_stdin(
            self.input.as_deref(),
        )?));
        let mut record = BedRegion::new();

        while input.next(&mut record)? {
            if record.start > expand_start {
                record.start -= expand_start;
            } else {
                record.start = 0;
            }
            if record.end < u64::MAX - expand_end {
                record.end += expand_end;
            } else {
                record.end = u64::MAX;
            }
            output.write_record(&record)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use clap::Parser;

    #[derive(Debug, Parser)]
    struct Cli {
        #[command(flatten)]
        extend_bed: ExtendBed,
    }

    #[test]
    fn test_expand() -> anyhow::Result<()> {
        let cli = Cli::parse_from(&[
            "command",
            "-o",
            "../target/expanded.bed",
            "--expand",
            "10",
            "testfiles/sample1.bed",
        ]);
        cli.extend_bed.run()?;

        let mut reader = BedReader::new(io::BufReader::new(autocompress::autodetect_open(
            "../target/expanded.bed",
        )?));

        let mut record = BedRegion::new();
        assert_eq!(reader.next(&mut record)?, true);
        assert_eq!(
            BedRegion {
                chromosome: b"1".to_vec(),
                start: 90,
                end: 210,
                columns: vec![]
            },
            record
        );

        Ok(())
    }

    #[test]
    fn test_expand2() -> anyhow::Result<()> {
        let cli = Cli::parse_from(&[
            "command",
            "-o",
            "../target/expanded2.bed",
            "--expand-start",
            "200",
            "--expand-end",
            "100",
            "testfiles/sample1.bed",
        ]);
        cli.extend_bed.run()?;

        let mut reader = BedReader::new(io::BufReader::new(autocompress::autodetect_open(
            "../target/expanded2.bed",
        )?));

        let mut record = BedRegion::new();
        assert_eq!(reader.next(&mut record)?, true);
        assert_eq!(
            BedRegion {
                chromosome: b"1".to_vec(),
                start: 0,
                end: 300,
                columns: vec![]
            },
            record
        );

        assert_eq!(reader.next(&mut record)?, true);
        assert_eq!(
            BedRegion {
                chromosome: b"1".to_vec(),
                start: 100,
                end: 500,
                columns: vec![]
            },
            record
        );

        Ok(())
    }
}
