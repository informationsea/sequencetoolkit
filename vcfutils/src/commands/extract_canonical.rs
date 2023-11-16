use clap::Args;
use std::io;
use std::io::Write;

#[derive(Debug, Args)]
#[command(
    about = "Extract canonical transcript from snpeff log",
    version,
    author
)]
pub struct ExtractCanonical {
    #[arg(
        help = "snpEff standard error log file.",
        long_help = r#"snpEff verbose standard error log file.

Please create snpEff standard error file with a command shown in below.
$ java -jar snpEff.jar ann -v -canon -noStats -noLog hg19 < /dev/null 2> snpeff.log
"#
    )]
    input: Option<String>,
    #[arg(short, long, help = "List of canonical transcripts")]
    output: Option<String>,
}

impl ExtractCanonical {
    //     fn command_name(&self) -> &'static str {
    //         "extract-canonical"
    //     }
    //     fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //         app.about("Extract canonical transcript from snpeff log")
    //             .arg(
    //                 Arg::with_name("input")
    //                     .index(1)
    //                     .takes_value(true)
    //                     .help("snpEff standard error log file.")
    //                     .long_help(
    //                         r#"snpEff verbose standard error log file.

    // Please create snpEff standard error file with a command shown in below.
    // $ java -jar snpEff.jar ann -v -canon -noStats -noLog hg19 < /dev/null 2> snpeff.log
    // "#,
    //                     ),
    //             )
    //             .arg(
    //                 Arg::with_name("output")
    //                     .short("o")
    //                     .long("output")
    //                     .takes_value(true)
    //                     .help("List of canonical transcripts"),
    //             )
    //     }

    pub fn run(&self) -> anyhow::Result<()> {
        let reader = io::BufReader::new(autocompress::autodetect_open_or_stdin(
            self.input.as_deref(),
        )?);
        let mut writer = autocompress::autodetect_create_or_stdout_prefer_bgzip(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?;
        let canonical_list = extract_canonical(reader)?;
        for one in canonical_list {
            writeln!(writer, "{}", one)?;
        }

        Ok(())
    }
}

fn extract_canonical(mut reader: impl io::BufRead) -> anyhow::Result<Vec<String>> {
    let mut result = Vec::new();

    let mut line = String::new();
    let mut in_list = false;
    while reader.read_line(&mut line)? > 0 {
        if line.starts_with("\t\tgeneName\tgeneId\ttranscriptId\tcdsLength") {
            in_list = true;
        } else if in_list && line.starts_with("\t\t") {
            let elements: Vec<_> = line.trim().split('\t').collect();
            if let Some(transcript) = elements.get(2) {
                result.push((*transcript).to_string());
            }
        }
        line.clear();
    }

    Ok(result)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_extract_canonical() -> anyhow::Result<()> {
        let snpeff_log = include_bytes!("../../testfiles/snpeff-log.txt");
        let canonical_list = extract_canonical(&snpeff_log[..])?;
        assert!(canonical_list.contains(&"ENST00000291182.9_4".to_string()));
        assert!(canonical_list.contains(&"ENST00000533876.1_1".to_string()));
        Ok(())
    }
}
