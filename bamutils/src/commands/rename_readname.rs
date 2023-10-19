use clap::Args;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::str;

#[derive(Debug, Args)]
#[command(about = "Rename read names", version, author)]
pub struct RenameReadname {
    #[arg(help = "Input BAM/CRAM file")]
    bam: String,
    #[arg(short = 'T', long, help = "Reference FASTA")]
    reference: Option<String>,
    #[arg(short, long, help = "Output BAM/CRAM/SAM file")]
    output: String,
    #[arg(short, long, help = "# of threads to write output file")]
    threads: Option<usize>,
}

impl RenameReadname {
    // fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //     app.about("Rename read name")
    //         .arg(
    //             Arg::with_name("bam")
    //                 .index(1)
    //                 .takes_value(true)
    //                 .required(true)
    //                 .help("Input BAM/CRAM file"),
    //         )
    //         .arg(
    //             Arg::with_name("reference")
    //                 .short("T")
    //                 .long("reference")
    //                 .takes_value(true)
    //                 .help("Reference FASTA"),
    //         )
    //         .arg(
    //             Arg::with_name("output")
    //                 .short("o")
    //                 .long("output")
    //                 .required(true)
    //                 .takes_value(true),
    //         )
    // }

    pub fn run(&self) -> anyhow::Result<()> {
        run(
            &self.bam,
            &self.output,
            self.reference.as_deref(),
            self.threads,
        )?;
        Ok(())
    }
}

fn run(
    bam_path: &str,
    output_path: &str,
    reference: Option<&str>,
    threads: Option<usize>,
) -> anyhow::Result<()> {
    let mut reader = bam::Reader::from_path(bam_path)?;
    if let Some(reference) = reference {
        reader.set_reference(reference)?;
    }

    let output_format = if output_path.ends_with(".cram") {
        bam::Format::Cram
    } else if output_path.ends_with(".bam") {
        bam::Format::Bam
    } else if output_path.ends_with(".sam") {
        bam::Format::Sam
    } else {
        return Err(anyhow::anyhow!("Failed to detect output format: {}", output_path).into());
    };
    let mut output_header = bam::Header::from_template(reader.header());
    output_header.push_comment(b"rename read names with sequence toolkit");
    let mut writer = bam::Writer::from_path(output_path, &output_header, output_format)?;
    if let Some(threads) = threads {
        writer.set_threads(threads)?;
    }
    if let Some(reference) = reference {
        writer.set_reference(reference)?;
    }

    let mut rename_readname: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut next_readname: u64 = 1;

    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        if let Some(d) = rename_readname.get(record.qname()) {
            record.set_qname(format!("R{}", d).as_bytes());
        } else {
            rename_readname.insert(record.qname().to_vec(), next_readname);
            record.set_qname(format!("R{}", next_readname).as_bytes());
            next_readname += 1;
        }
        writer.write(&record)?;
    }

    Ok(())
}
