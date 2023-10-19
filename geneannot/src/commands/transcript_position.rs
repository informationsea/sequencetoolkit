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
pub struct TranscriptPosition {
    #[arg(
        help = "geneannot database (INPUT / gzip BINCODE)",
        short = 'd',
        long = "database"
    )]
    db: String,
    #[arg(help = "Genome position", short = 'p', long = "position")]
    position: String,
    #[arg(help = "Chromosome name", short = 'c', long = "chromosome")]
    chromosome: String,
}

impl TranscriptPosition {
    // fn command_name(&self) -> &'static str {
    //     "transcript-position"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("Convert CDS/transcript position into genome position")
    //         .arg(
    //             Arg::with_name("db")
    //                 .help("geneannot database (INPUT / gzip BINCODE)")
    //                 .short("d")
    //                 .long("database")
    //                 .required(true)
    //                 .takes_value(true),
    //         )
    //         .arg(
    //             Arg::with_name("position")
    //                 .help("Genome position")
    //                 .short("p")
    //                 .long("position")
    //                 .required(true)
    //                 .takes_value(true),
    //         )
    //         .arg(
    //             Arg::with_name("chromosome")
    //                 .help("Chromosome name")
    //                 .short("c")
    //                 .long("chromosome")
    //                 .required(true)
    //                 .takes_value(true),
    //         )
    // }
    pub fn run(&self) -> anyhow::Result<()> {
        Ok(search_transcript_position(
            &self.db,
            &self.chromosome,
            &self.position,
        )?)
    }
}

fn search_transcript_position(
    db: &str,
    chromosome: &str,
    position: &str,
) -> Result<(), GeneAnnotError> {
    let db_reader = BufReader::new(MultiGzDecoder::new(File::open(db)?));
    let db: crate::annotator::models::GeneAnnotations = bincode::deserialize_from(db_reader)?;
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
        Err(GeneAnnotError::OtherError("Unknown chromosome name").into())
    }
}
