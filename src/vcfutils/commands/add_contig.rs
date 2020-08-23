use super::Command;
use crate::vcfutils::logic::add_contig;
use clap::{App, Arg, ArgMatches};
use std::io::BufReader;

pub struct AddContig;

impl Command for AddContig {
    fn command_name(&self) -> &'static str {
        "add-contig"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Scan VCF file and add contig header")
            .arg(
                Arg::with_name("input")
                    .index(1)
                    .takes_value(true)
                    .required(true)
                    .help("Input VCF file"),
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("output")
                    .takes_value(true)
                    .help("Output VCF"),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        let vcf_reader = autocompress::open(matches.value_of("input").unwrap())?;
        let mut vcf_writer = autocompress::create_or_stdout(matches.value_of("output"))?;

        let contigs = add_contig::scan_contig(&mut BufReader::new(vcf_reader))?;
        let vcf_reader = autocompress::open(matches.value_of("input").unwrap())?;
        add_contig::add_contig(&mut BufReader::new(vcf_reader), &mut vcf_writer, &contigs)?;

        Ok(())
    }
}
