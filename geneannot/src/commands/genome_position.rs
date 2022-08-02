use super::Command;
use crate::annotator;
use crate::annotator::hgvs::position::{parse_hgvs_position, ParsedPosition};
use crate::annotator::models::TranscriptTrait;
use crate::GeneAnnotError;
use clap::{App, Arg, ArgMatches};
use flate2::read::MultiGzDecoder;
use log::info;
use std::fs::File;
use std::io::BufReader;

pub struct GenomePosition;

impl Command for GenomePosition {
    fn command_name(&self) -> &'static str {
        "genome-position"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Convert genome position into CDS/transcript position")
            .arg(
                Arg::with_name("db")
                    .help("geneannot database (INPUT / gzip BINCODE)")
                    .short("d")
                    .long("database")
                    .required(true)
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("position")
                    .help("CDS or transcript position in HGVS Sequence Variant Nomenclature")
                    .short("p")
                    .long("position")
                    .required(true)
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("transcript-name")
                    .help("Transcript name")
                    .short("t")
                    .long("transcript")
                    .required(true)
                    .takes_value(true),
            )
    }
    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        Ok(search_genome_position(
            matches.value_of("db").unwrap(),
            matches.value_of("transcript-name").unwrap(),
            matches.value_of("position").unwrap(),
        )?)
    }
}

fn search_genome_position(
    db: &str,
    transcript_name: &str,
    position: &str,
) -> Result<(), GeneAnnotError> {
    let db_reader = BufReader::new(MultiGzDecoder::new(File::open(db)?));
    let db: annotator::models::GeneAnnotations = bincode::deserialize_from(db_reader)?;
    info!("database loaded");
    let parsed_position = parse_hgvs_position(position)?;

    if let Some((_, transcript)) = db.transcript(transcript_name) {
        match parsed_position {
            ParsedPosition::GenomePosition(_) => {
                return Err(GeneAnnotError::HgvsPositionParseError);
            }
            ParsedPosition::CdsPosition(x) => {
                if let Some(p) = transcript.genome_position_from_cds(x) {
                    println!(
                        "{}:g.{}",
                        db.genome().chromosomes()[transcript.chromosome_index()].name,
                        p + 1
                    );
                } else {
                    return Err(GeneAnnotError::OtherError("No CDS"));
                }
            }
            ParsedPosition::TranscriptPosition(x) => {
                let p = transcript.genome_position(x);
                println!(
                    "{}:g.{}",
                    db.genome().chromosomes()[transcript.chromosome_index()].name,
                    p + 1
                );
            }
        }
        Ok(())
    } else {
        Err(GeneAnnotError::OtherError("Transcript is not found"))
    }
}
