use crate::utils::fastq::*;
use crate::utils::largereorder::LargeReorder;
use crate::utils::{DNBSeqFASTQInfo, DigestWriter};
use anyhow::Context;
use autocompress::io::RayonReader;
use clap::{Args, ValueEnum};
use std::sync::mpsc::channel;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Read {
    #[value(name = "1")]
    One,
    #[value(name = "2")]
    Two,
}

#[derive(Debug, Clone, Args)]
#[command(about = "Reorder DNBSeq FASTQ", version, author)]
pub struct DNBSeqFastqReorder {
    #[arg(help = "Input FASTQ")]
    input: String,
    #[arg(short, long, help = "Output FASTQ")]
    output: String,
    #[arg(help = "Original FASTQ information", short = 'i', long = "fastq-info")]
    fastq_info: String,
    #[arg(help = "Read number", short = 'r', long = "read")]
    read: Option<Read>,
    #[arg(
        help = "# of threads",
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

impl DNBSeqFastqReorder {
    pub fn run(&self) -> Result<(), anyhow::Error> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.thread)
            .build_global()
            .context("Failed to set # of threads")?;

        let mut input_fastq = RayonReader::new(
            autocompress::autodetect_open(&self.input).context("Failed to open input FASTQ")?,
        );

        let mut output = DigestWriter::new(
            autocompress::autodetect_parallel_create(
                &self.output,
                autocompress::CompressionLevel::Default,
            )?,
            md5::Md5::default(),
        );

        let fastq_info = DNBSeqFASTQInfo::load(RayonReader::new(autocompress::autodetect_open(
            self.fastq_info.as_str(),
        )?))?;

        let (sender, receiver) = channel();

        let mut fastq_entries = LargeReorder::with_temp_dir(
            self.on_memory_count,
            self.temporary_directory.as_deref(),
            Some("dnbseqfastq."),
        );
        let read = self.read;

        rayon::spawn(move || {
            let f = move || -> anyhow::Result<()> {
                eprintln!("Loading FASTQ...");
                while let Some(entry) = GenericFastqEntry::read(&mut input_fastq)? {
                    fastq_entries.add(DNBSeqFastqEntry::from_fastq_entry(entry, &fastq_info)?)?;
                }

                eprintln!("Writing FASTQ...");
                for one in fastq_entries.into_iter()? {
                    let one = one?;
                    one.write(
                        &mut output,
                        match read {
                            Some(Read::One) => Some("1"),
                            Some(Read::Two) => Some("2"),
                            None => None,
                        },
                        &fastq_info,
                    )?;
                }

                let md5result = format!("{:x}", output.finalize_reset());
                eprintln!("MD5: {}", md5result);

                if let Some(read) = read {
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
                Ok(())
            };
            sender.send(f()).expect("Failed to send final result");
        });

        receiver.recv().expect("Failed to receive final result")
    }
}
