use crate::annotator;
use crate::annotator::hgvs::position::{parse_hgvs_position, ParsedPosition};
use crate::annotator::models::TranscriptTrait;
use crate::GeneAnnotError;
use clap::Args;
use flate2::read::MultiGzDecoder;
use log::info;
use std::fs::File;
use std::io::BufReader;

#[derive(Debug, Args)]
#[command(
    about = "Convert CDS/transcript position into genome position",
    version,
    author
)]
pub struct GenomePosition {
    #[arg(
        help = "geneannot database (INPUT / gzip BINCODE)",
        short = 'd',
        long = "database"
    )]
    db: String,
    #[arg(
        help = "CDS or transcript position in HGVS Sequence Variant Nomenclature",
        short = 'p',
        long = "position"
    )]
    position: String,
    #[arg(help = "Transcript name", short = 't', long = "transcript-name")]
    transcript_name: String,
}

impl GenomePosition {
    pub fn run(&self) -> anyhow::Result<()> {
        Ok(search_genome_position(
            &self.db,
            &self.transcript_name,
            &self.position,
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
