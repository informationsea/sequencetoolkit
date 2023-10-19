use bio::io::fasta::IndexedReader;
use clap::Args;
use flate2::{write::GzEncoder, Compression};
use std::fs::File;
use std::path::Path;

#[derive(Debug, Args)]
#[command(
    about = "Create geneannot database from refGene format database",
    version,
    author
)]
pub struct CreateDb {
    #[arg(
        help = "refGene format database (INPUT / plain text or gzip)",
        long_help = r#"refGene format database (INPUT / plain text or gzip)
example file URL:
- http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
- http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV33.txt.gz"#
    )]
    db: String,
    #[arg(
        short = 'o',
        long,
        help = "database output path (OUTPUT / gzip BINCODE)"
    )]
    output: String,
    #[arg(
        short = 'f',
        long,
        help = "Reference FASTA file (INPUT / indexed FASTA)"
    )]
    fasta: String,
    #[arg(short = 'n', long, help = "Reference name")]
    reference_name: Option<String>,
}

impl CreateDb {
    //     fn command_name(&self) -> &'static str {
    //         "create-db"
    //     }
    //     fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //         app.about("Create gene database from refGene format file")
    //             .arg(
    //                 Arg::with_name("db")
    //                     .help("refGene format database (INPUT / plain text or gzip)")
    //                     .short("r")
    //                     .long("refgene")
    //                     .required(true)
    //                     .takes_value(true)
    //                     .long_help(
    //                         r#"refGene format database (INPUT / plain text or gzip)
    // example file URL:
    //   - http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
    //   - http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV33.txt.gz
    // "#,
    //                     ),
    //             )
    //             .arg(
    //                 Arg::with_name("output")
    //                     .help("database output path (OUTPUT / gzip BINCODE)")
    //                     .short("o")
    //                     .long("output")
    //                     .required(true)
    //                     .takes_value(true),
    //             )
    //             .arg(
    //                 Arg::with_name("fasta")
    //                     .help("Reference FASTA file (INPUT / indexed FASTA)")
    //                     .short("f")
    //                     .long("reference")
    //                     .required(true)
    //                     .takes_value(true),
    //             )
    //             .arg(
    //                 Arg::with_name("reference-name")
    //                     .help("Reference name")
    //                     .short("n")
    //                     .long("reference-name")
    //                     .takes_value(true),
    //             )
    //     }

    pub fn run(&self) -> anyhow::Result<()> {
        Ok(create_db(
            &self.db,
            &self.output,
            &self.fasta,
            self.reference_name.as_deref(),
        )?)
    }
}

fn create_db(
    ref_gene: &str,
    output: &str,
    fasta: &str,
    reference_name: Option<&str>,
) -> Result<(), crate::GeneAnnotError> {
    let db_reader = autocompress::open(ref_gene)?;
    let output_writer = GzEncoder::new(File::create(output)?, Compression::default());
    let fasta_reader = IndexedReader::from_file(&fasta)?;
    let reference_name = reference_name.unwrap_or_else(|| {
        Path::new(fasta)
            .file_name()
            .map(|x| x.to_str().unwrap())
            .unwrap_or("reference")
    });
    let genome = crate::annotator::models::load_fasta(reference_name, &fasta_reader.index);
    let gene_annotation = crate::annotator::models::load_refgene(genome, db_reader)?;
    //serde_json::to_writer(output_writer, &gene_annotation)?;
    bincode::serialize_into(output_writer, &gene_annotation)?;

    Ok(())
}
