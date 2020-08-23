use super::Command;
use crate::geneannot::annotator;
use crate::geneannot::annotator::models::TranscriptTrait;
use crate::geneannot::{GeneAnnotError, GeneAnnotErrorKind};
use clap::{App, Arg, ArgMatches};
use flate2::read::MultiGzDecoder;
use log::info;
use std::fs::File;
use std::io::BufReader;

pub struct TranscriptPosition;

impl Command for TranscriptPosition {
    fn command_name(&self) -> &'static str {
        "transcript-position"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Convert CDS/transcript position into genome position")
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
                    .help("Genome position")
                    .short("p")
                    .long("position")
                    .required(true)
                    .takes_value(true),
            )
            .arg(
                Arg::with_name("chromosome")
                    .help("Chromosome name")
                    .short("c")
                    .long("chromosome")
                    .required(true)
                    .takes_value(true),
            )
    }
    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::SequenceToolkitError> {
        Ok(search_transcript_position(
            matches.value_of("db").unwrap(),
            matches.value_of("chromosome").unwrap(),
            matches.value_of("position").unwrap(),
        )?)
    }
}

fn search_transcript_position(
    db: &str,
    chromosome: &str,
    position: &str,
) -> Result<(), GeneAnnotError> {
    let db_reader = BufReader::new(MultiGzDecoder::new(File::open(db)?));
    let db: annotator::models::GeneAnnotations = bincode::deserialize_from(db_reader)?;
    info!("database loaded");
    let position = position.parse::<u64>()? - 1;
    if let Some(chromosome_index) = db.genome().chromosome_index(chromosome) {
        for one in db
            .interval_tree(chromosome_index)
            .unwrap()
            .find(position..(position + 1))
        {
            let gene = &db.genes()[one.data().0];
            let transcript = &gene.transcripts()[one.data().1];
            print!("{}({}):", transcript.id(), gene.id());
            if let Some(cds_position) = transcript.cds_position(position) {
                println!("{}", cds_position);
            } else {
                println!("{}", transcript.transcript_position(position));
            }
        }
        Ok(())
    } else {
        Err(GeneAnnotErrorKind::OtherError("Unknown chromosome name").into())
    }
}
