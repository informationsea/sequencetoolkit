use crate::utils::fastq::*;
use crate::utils::largereorder::LargeReorder;
use crate::utils::{DigestWriter, IlluminaFASTQInfo};
use anyhow::Context;
use std::convert::TryFrom;
use std::io;
use std::str;

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum Read {
    #[value(name = "1")]
    One,
    #[value(name = "2")]
    Two,
}

impl Read {
    pub fn to_str(self) -> &'static str {
        match self {
            Read::One => "1",
            Read::Two => "2",
        }
    }
}

#[derive(Debug, Clone, clap::Args)]
#[command(about = "Reorder Illumina FASTQ", version, author)]
pub struct IlluminaFastqReorder {
    #[arg(help = "Input FASTQ")]
    input: String,
    #[arg(help = "Output", short = 'o', long = "output")]
    output: String,
    #[arg(help = "Original FASTQ information", short = 'i', long = "fastq-info")]
    fastq_info: Option<String>,
    #[arg(help = "Read number", short = 'r', long = "read")]
    read: Option<Read>,
    #[arg(
        help = "# of write threads",
        short = 't',
        long = "thread",
        default_value = "1"
    )]
    thread: usize,
    #[arg(help = "Temporary files directory", short = 't', long = "temporary")]
    temporary_directory: Option<String>,
    #[arg(
        help = "# of entries on the memory",
        short = 'm',
        long = "mem",
        default_value = "1000000"
    )]
    on_memory_count: usize,
}

impl IlluminaFastqReorder {
    pub fn run(&self) -> Result<(), anyhow::Error> {
        let iothread = autocompress::iothread::IoThread::new(self.thread);

        let mut input_fastq = iothread
            .open(&self.input)
            .context("Cannot open input FASTQ")?;

        let mut output = DigestWriter::new(
            iothread
                .create(&self.output, autocompress::CompressionLevel::Default)
                .context("Cannot open output FASTQ")?,
            md5::Md5::default(),
        );

        let fastq_info = if let Some(fastq_info_path) = self.fastq_info.as_deref() {
            Some(IlluminaFASTQInfo::load(io::BufReader::new(
                autocompress::open(fastq_info_path).context("Cannot read FASTQ Info")?,
            ))?)
        } else {
            None
        };

        let mut fastq_entries = LargeReorder::with_temp_dir(
            self.on_memory_count,
            self.temporary_directory.as_deref(),
            Some("illuminafastq."),
        );

        eprintln!("Loading FASTQ...");
        while let Some(entry) = GenericFastqEntry::read(&mut input_fastq)? {
            fastq_entries.add(IlluminaFastqEntry::try_from(entry)?)?;
        }

        eprintln!("Writing FASTQ...");
        for one in fastq_entries.into_iter()? {
            let one = one?;
            one.write(
                &mut output,
                self.read.map(|x| x.to_str()),
                fastq_info.as_ref(),
            )?;
        }

        let md5result = format!("{:x}", output.finalize_reset());
        eprintln!("MD5: {}", md5result);

        if let Some(fastq_info) = fastq_info {
            if let Some(read) = self.read {
                match read {
                    Read::One => {
                        if fastq_info.fastq1_md5 == md5result {
                            eprintln!("MD5 OK");
                        } else {
                            eprintln!("MD5 NG. Expected MD5 is {}", fastq_info.fastq1_md5);
                        }
                    }
                    Read::Two => {
                        if fastq_info.fastq2_md5 == md5result {
                            eprintln!("MD5 OK");
                        } else {
                            eprintln!("MD5 NG. Expected MD5 is {}", fastq_info.fastq2_md5);
                        }
                    }
                }
            }
        }

        Ok(())
    }
}
