use super::Command;
use crate::vcfutils::logic::replace_contig::replace_contig;
use crate::vcfutils::utils;
use clap::{App, Arg, ArgMatches};
use std::collections::HashMap;
use std::io::BufReader;
use vcf::U8Vec;

pub struct ReplaceContig;

impl Command for ReplaceContig {
    fn command_name(&self) -> &'static str {
        "replace-contig"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Replace contig")
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
                    .help("Output file"),
            )
            .arg(
                Arg::with_name("contig-mapping")
                    .short("m")
                    .long("contig-mapping")
                    .takes_value(true)
                    .help("contig mapping file (csv/tsv)")
                    .required_unless_one(&[
                        "replace-refseq",
                        "add-chr-prefix",
                        "remove-chr-prefix",
                    ]),
            )
            .arg(
                Arg::with_name("replace-refseq")
                    .short("e")
                    .long("replace-refseq")
                    .help("Replace RefSeq contig name with chromosome number (for dbSNP 152/153)"),
            )
            .arg(
                Arg::with_name("add-chr-prefix")
                    .short("a")
                    .long("add-chr-prefix")
                    .help("add \"chr\" prefix to human chromosomes"),
            )
            .arg(
                Arg::with_name("remove-chr-prefix")
                    .short("r")
                    .long("remove-chr-prefix")
                    .help("remove \"chr\" prefix from human chromosomes"),
            )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        let mut reader = BufReader::new(autocompress::open_or_stdin(matches.value_of("input"))?);
        let mut writer = autocompress::create_or_stdout(matches.value_of("output"))?;
        let mapping = if matches.is_present("replace-refseq") {
            refseq_replace()
        } else if matches.is_present("add-chr-prefix") {
            add_chr()
        } else if matches.is_present("remove-chr-prefix") {
            remove_chr()
        } else {
            utils::load_mapping(utils::auto_csv_reader_from_path(
                matches.value_of("contig-mapping").unwrap(),
                false,
            )?)?
            .mapping
        };

        replace_contig(&mut reader, &mut writer, &mapping)?;
        Ok(())
    }
}

fn add_chr() -> HashMap<U8Vec, U8Vec> {
    [
        (b"1".to_vec(), b"chr1".to_vec()),
        (b"2".to_vec(), b"chr2".to_vec()),
        (b"3".to_vec(), b"chr3".to_vec()),
        (b"4".to_vec(), b"chr4".to_vec()),
        (b"5".to_vec(), b"chr5".to_vec()),
        (b"6".to_vec(), b"chr6".to_vec()),
        (b"7".to_vec(), b"chr7".to_vec()),
        (b"8".to_vec(), b"chr8".to_vec()),
        (b"9".to_vec(), b"chr9".to_vec()),
        (b"10".to_vec(), b"chr10".to_vec()),
        (b"11".to_vec(), b"chr11".to_vec()),
        (b"12".to_vec(), b"chr12".to_vec()),
        (b"13".to_vec(), b"chr13".to_vec()),
        (b"14".to_vec(), b"chr14".to_vec()),
        (b"15".to_vec(), b"chr15".to_vec()),
        (b"16".to_vec(), b"chr16".to_vec()),
        (b"17".to_vec(), b"chr17".to_vec()),
        (b"18".to_vec(), b"chr18".to_vec()),
        (b"19".to_vec(), b"chr19".to_vec()),
        (b"20".to_vec(), b"chr20".to_vec()),
        (b"21".to_vec(), b"chr21".to_vec()),
        (b"22".to_vec(), b"chr22".to_vec()),
        (b"X".to_vec(), b"chrX".to_vec()),
        (b"Y".to_vec(), b"chrY".to_vec()),
    ]
    .iter()
    .cloned()
    .collect()
}

fn remove_chr() -> HashMap<U8Vec, U8Vec> {
    [
        (b"chr1".to_vec(), b"1".to_vec()),
        (b"chr2".to_vec(), b"2".to_vec()),
        (b"chr3".to_vec(), b"3".to_vec()),
        (b"chr4".to_vec(), b"4".to_vec()),
        (b"chr5".to_vec(), b"5".to_vec()),
        (b"chr6".to_vec(), b"6".to_vec()),
        (b"chr7".to_vec(), b"7".to_vec()),
        (b"chr8".to_vec(), b"8".to_vec()),
        (b"chr9".to_vec(), b"9".to_vec()),
        (b"chr10".to_vec(), b"10".to_vec()),
        (b"chr11".to_vec(), b"11".to_vec()),
        (b"chr12".to_vec(), b"12".to_vec()),
        (b"chr13".to_vec(), b"13".to_vec()),
        (b"chr14".to_vec(), b"14".to_vec()),
        (b"chr15".to_vec(), b"15".to_vec()),
        (b"chr16".to_vec(), b"16".to_vec()),
        (b"chr17".to_vec(), b"17".to_vec()),
        (b"chr18".to_vec(), b"18".to_vec()),
        (b"chr19".to_vec(), b"19".to_vec()),
        (b"chr20".to_vec(), b"20".to_vec()),
        (b"chr21".to_vec(), b"21".to_vec()),
        (b"chr22".to_vec(), b"22".to_vec()),
        (b"chrX".to_vec(), b"X".to_vec()),
        (b"chrY".to_vec(), b"Y".to_vec()),
    ]
    .iter()
    .cloned()
    .collect()
}

fn refseq_replace() -> HashMap<U8Vec, U8Vec> {
    [
        (b"NC_000001.10".to_vec(), b"1".to_vec()),
        (b"NC_000002.11".to_vec(), b"2".to_vec()),
        (b"NC_000003.11".to_vec(), b"3".to_vec()),
        (b"NC_000004.11".to_vec(), b"4".to_vec()),
        (b"NC_000005.9".to_vec(), b"5".to_vec()),
        (b"NC_000006.11".to_vec(), b"6".to_vec()),
        (b"NC_000007.13".to_vec(), b"7".to_vec()),
        (b"NC_000008.10".to_vec(), b"8".to_vec()),
        (b"NC_000009.11".to_vec(), b"9".to_vec()),
        (b"NC_000010.10".to_vec(), b"10".to_vec()),
        (b"NC_000011.9".to_vec(), b"11".to_vec()),
        (b"NC_000012.11".to_vec(), b"12".to_vec()),
        (b"NC_000013.10".to_vec(), b"13".to_vec()),
        (b"NC_000014.8".to_vec(), b"14".to_vec()),
        (b"NC_000015.9".to_vec(), b"15".to_vec()),
        (b"NC_000016.9".to_vec(), b"16".to_vec()),
        (b"NC_000017.10".to_vec(), b"17".to_vec()),
        (b"NC_000018.9".to_vec(), b"18".to_vec()),
        (b"NC_000019.9".to_vec(), b"19".to_vec()),
        (b"NC_000020.10".to_vec(), b"20".to_vec()),
        (b"NC_000021.8".to_vec(), b"21".to_vec()),
        (b"NC_000022.10".to_vec(), b"22".to_vec()),
        (b"NC_000023.10".to_vec(), b"X".to_vec()),
        (b"NC_000024.9".to_vec(), b"Y".to_vec()),
        (b"NC_000001.11".to_vec(), b"1".to_vec()),
        (b"NC_000002.12".to_vec(), b"2".to_vec()),
        (b"NC_000003.12".to_vec(), b"3".to_vec()),
        (b"NC_000004.12".to_vec(), b"4".to_vec()),
        (b"NC_000005.10".to_vec(), b"5".to_vec()),
        (b"NC_000006.12".to_vec(), b"6".to_vec()),
        (b"NC_000007.14".to_vec(), b"7".to_vec()),
        (b"NC_000008.11".to_vec(), b"8".to_vec()),
        (b"NC_000009.12".to_vec(), b"9".to_vec()),
        (b"NC_000010.11".to_vec(), b"10".to_vec()),
        (b"NC_000011.10".to_vec(), b"11".to_vec()),
        (b"NC_000012.12".to_vec(), b"12".to_vec()),
        (b"NC_000013.11".to_vec(), b"13".to_vec()),
        (b"NC_000014.9".to_vec(), b"14".to_vec()),
        (b"NC_000015.10".to_vec(), b"15".to_vec()),
        (b"NC_000016.10".to_vec(), b"16".to_vec()),
        (b"NC_000017.11".to_vec(), b"17".to_vec()),
        (b"NC_000018.10".to_vec(), b"18".to_vec()),
        (b"NC_000019.10".to_vec(), b"19".to_vec()),
        (b"NC_000020.11".to_vec(), b"20".to_vec()),
        (b"NC_000021.9".to_vec(), b"21".to_vec()),
        (b"NC_000022.11".to_vec(), b"22".to_vec()),
        (b"NC_000023.11".to_vec(), b"X".to_vec()),
        (b"NC_000024.10".to_vec(), b"Y".to_vec()),
    ]
    .iter()
    .cloned()
    .collect()
}
