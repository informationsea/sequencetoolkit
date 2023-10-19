use super::vcf2csv::create_config;
use crate::commands::vcf2csv::TableConfig;
use crate::logic::generate_sql::generate_sql;
use crate::utils;
use clap::Args;
use std::str;

#[derive(Debug, Args)]
#[command(about = "Generate SQL from VCF", version, author)]
pub struct GenerateSql {
    #[arg(help = "Input VCF file")]
    input: Option<String>,
    #[arg(short, long, help = "Output SQL")]
    output: Option<String>,
    #[arg(short, long, help = "INFO tags to include")]
    info: Option<Vec<String>>,
    #[arg(short, long, help = "FORMAT tags to include")]
    format: Option<Vec<String>>,
    #[arg(short = 'n', long, help = "Table name")]
    table_name: String,
}

impl TableConfig for GenerateSql {
    fn canonical_list(&self) -> Option<&str> {
        None
    }
    fn input(&self) -> &[String] {
        &[] // TODO: FIX HERE
    }
    fn replace_sample_name(&self) -> Option<&[String]> {
        None
    }
    fn group_names(&self) -> Option<&[String]> {
        None
    }
    fn split_multi_allelic(&self) -> bool {
        false
    }
    fn decode_genotype(&self) -> bool {
        false
    }
    fn info_list(&self) -> Option<&[String]> {
        self.info.as_deref()
    }
    fn format_list(&self) -> Option<&[String]> {
        self.format.as_deref()
    }
}

impl GenerateSql {
    // fn command_name(&self) -> &'static str {
    //     "generate-sql"
    // }
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("List up sample names")
    //         .arg(
    //             Arg::with_name("input")
    //                 .index(1)
    //                 .takes_value(true)
    //                 .help("Input VCF file"),
    //         )
    //         .arg(
    //             Arg::with_name("output")
    //                 .short("o")
    //                 .long("output")
    //                 .takes_value(true)
    //                 .help("Output SQL"),
    //         )
    //         .arg(
    //             Arg::with_name("info")
    //                 .short("i")
    //                 .long("info")
    //                 .help("INFO tags to include")
    //                 .takes_value(true)
    //                 .multiple(true),
    //         )
    //         .arg(
    //             Arg::with_name("format")
    //                 .short("f")
    //                 .long("format")
    //                 .help("FORMAT tags to include")
    //                 .takes_value(true)
    //                 .multiple(true),
    //         )
    //         .arg(
    //             Arg::with_name("table-name")
    //                 .short("n")
    //                 .long("table-name")
    //                 .help("Table name")
    //                 .takes_value(true)
    //                 .required(true),
    //         )
    // }

    pub fn run(&self) -> anyhow::Result<()> {
        let vcf_reader = utils::open_vcf_from_path(self.input.as_deref())?;
        let mut output = autocompress::create_or_stdout(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?;
        let mut config = create_config(&vcf_reader.header(), self)?;
        config.split_multi_allelic = true;
        config.decoded_genotype = false;
        let generated_sql = generate_sql(&vcf_reader.header(), &config, &self.table_name)?;
        output.write_all(generated_sql.as_bytes())?;
        Ok(())
    }
}
