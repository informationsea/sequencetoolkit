use super::super::{BedReader, BedRegion, BedWriter};
use crate::Command;
use clap::{App, Arg, ArgMatches};
use failure::ResultExt;
use sequencetoolkit_common::SequenceToolkitErrorKind;
use std::io;
use std::str;
use std::u64;

pub struct ExtendBed;

impl Command for ExtendBed {
    fn command_name(&self) -> &'static str {
        "expand"
    }

    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Extend BED regions")
            .arg(
                Arg::with_name("input")
                    .index(1)
                    .takes_value(true)
                    .help("Input BED files"),
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("output")
                    .takes_value(true)
                    .help("Output Text"),
            )
            .arg(
                Arg::with_name("expand")
                    .short("e")
                    .long("expand")
                    .takes_value(true)
                    .help("Expand length")
                    .required_unless("expand-start"),
            )
            .arg(
                Arg::with_name("expand-start")
                    .short("s")
                    .long("expand-start")
                    .takes_value(true)
                    .help("Expand start position")
                    .conflicts_with("expand")
                    .requires("expand-end"),
            )
            .arg(
                Arg::with_name("expand-end")
                    .short("n")
                    .long("expand-end")
                    .takes_value(true)
                    .help("Expand end position")
                    .conflicts_with("expand")
                    .requires("expand-start"),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        let expand_start: u64 = matches
            .value_of("expand-start")
            .or_else(|| matches.value_of("expand"))
            .unwrap()
            .parse::<u64>()
            .context(SequenceToolkitErrorKind::OtherError(
                "Cannot parse expand length",
            ))?;
        let expand_end: u64 = matches
            .value_of("expand-end")
            .or_else(|| matches.value_of("expand"))
            .unwrap()
            .parse::<u64>()
            .context(SequenceToolkitErrorKind::OtherError(
                "Cannot parse expand length",
            ))?;

        let mut output = BedWriter::new(io::BufWriter::new(autocompress::create_or_stdout(
            matches.value_of("output"),
        )?));
        let mut input = BedReader::new(io::BufReader::new(autocompress::open_or_stdin(
            matches.value_of("input"),
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

    #[test]
    fn test_expand() -> Result<(), crate::SequenceToolkitError> {
        let app = ExtendBed {};
        let command = app.config_subcommand(App::new("test"));
        let matches = command.get_matches_from(&[
            "command",
            "-o",
            "../target/expanded.bed",
            "--expand",
            "10",
            "testfiles/sample1.bed",
        ]);
        app.run(&matches)?;

        let mut reader = BedReader::new(io::BufReader::new(autocompress::open(
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
    fn test_expand2() -> Result<(), crate::SequenceToolkitError> {
        let app = ExtendBed {};
        let command = app.config_subcommand(App::new("test"));
        let matches = command.get_matches_from(&[
            "command",
            "-o",
            "../target/expanded2.bed",
            "--expand-start",
            "200",
            "--expand-end",
            "100",
            "testfiles/sample1.bed",
        ]);
        app.run(&matches)?;

        let mut reader = BedReader::new(io::BufReader::new(autocompress::open(
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
