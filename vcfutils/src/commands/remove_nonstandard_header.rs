use super::Command;
use clap::{App, Arg, ArgMatches};
use std::io::{self, BufRead};
use std::str;

pub struct RemoveNonStandardHeader;

impl Command for RemoveNonStandardHeader {
    fn command_name(&self) -> &'static str {
        "remove-non-standard-header"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Remove non standard header and reorder header")
            .arg(
                Arg::with_name("input")
                    .index(1)
                    .takes_value(true)
                    .help("Input VCF file"),
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("output")
                    .takes_value(true)
                    .help("Output CSV"),
            )
            .arg(
                Arg::with_name("additional")
                    .short("a")
                    .long("additional-header")
                    .takes_value(true)
                    .multiple(true)
                    .help("Add additional header line"),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        let mut reader =
            io::BufReader::new(autocompress::open_or_stdin(matches.value_of("input"))?);
        let mut writer = autocompress::create_or_stdout(
            matches.value_of("output"),
            autocompress::CompressionLevel::Default,
        )?;

        let mut header_file_format: Vec<Vec<u8>> = Vec::new();
        let mut header_info: Vec<Vec<u8>> = Vec::new();
        let mut header_format: Vec<Vec<u8>> = Vec::new();
        let mut header_alt: Vec<Vec<u8>> = Vec::new();
        let mut header_filter: Vec<Vec<u8>> = Vec::new();
        let mut header_contig: Vec<Vec<u8>> = Vec::new();

        let mut line: Vec<u8> = Vec::new();
        let mut index = 0;
        while reader.read_until(b'\n', &mut line)? > 0 {
            index += 1;
            if line.starts_with(b"##") {
                let header_line = vcf::VCFHeaderLine::from_bytes(&line, index)?;
                match header_line.contents() {
                    vcf::VCFHeaderContent::ALT { .. } => {
                        header_alt.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::INFO { .. } => {
                        header_info.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::FILTER { .. } => {
                        header_filter.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::FileFormat { .. } => {
                        header_file_format.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::Contig { .. } => {
                        header_contig.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::FORMAT { .. } => {
                        header_format.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::Other { .. } => (),
                }
            } else if line.starts_with(b"#") {
                if let Some(additional) = matches.values_of("additional") {
                    for one in header_file_format {
                        writer.write_all(&one)?;
                    }
                    for one in additional {
                        writer.write_all(one.as_bytes())?;
                        writer.write_all(b"\n")?;
                    }
                    for one_type in [
                        header_alt,
                        header_filter,
                        header_info,
                        header_format,
                        header_contig,
                    ]
                    .iter()
                    {
                        for one_line in one_type {
                            writer.write_all(one_line)?;
                        }
                    }
                }
                writer.write_all(&line)?;
                break;
            } else {
                writer.write_all(&line)?;
                break;
            }
            line.clear();
        }

        io::copy(&mut reader, &mut writer)?;

        Ok(())
    }
}
