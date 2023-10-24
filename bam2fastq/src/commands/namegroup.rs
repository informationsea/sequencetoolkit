use clap::{Args, ValueEnum};
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Format {
    Auto,
    SAM,
    BAM,
    CRAM,
}

#[derive(Debug, Clone, Args)]
#[command(about = "Create name grouped SAM/BAM/CRAM", version, author)]
pub struct NameGroup {
    #[arg(help = "Input SAM/BAM/CRAM")]
    input: Option<String>,
    #[arg(help = "Output SAM/BAM/CRAM", short = 'o', long = "output")]
    output: Option<String>,
    #[arg(
        help = "Output format (Supplementary and secondary mapped reads will be removed)",
        short = 'O',
        long = "output-format",
        default_value = "Auto"
    )]
    output_format: Format,
    #[arg(long = "read-threads")]
    read_threads: Option<usize>,
    #[arg(long = "write-threads")]
    write_threads: Option<usize>,
    #[arg(long = "compression-level")]
    compression_level: Option<u8>,
    #[arg(long = "reference", short = 'T')]
    reference: Option<String>,
}

impl NameGroup {
    pub fn run(&self) -> Result<(), anyhow::Error> {
        let mut input_bam = if let Some(bam_path) = self.input.as_ref() {
            if bam_path.starts_with("http://")
                || bam_path.starts_with("https://")
                || bam_path.starts_with("s3://")
            {
                bam::Reader::from_url(&url::Url::parse(bam_path)?)
            } else {
                bam::Reader::from_path(bam_path)
            }
        } else {
            bam::Reader::from_stdin()
        }?;

        if let Some(reference) = self.reference.as_deref() {
            input_bam.set_reference(reference)?;
        }

        if let Some(read_threads) = self.read_threads {
            input_bam.set_threads(read_threads)?;
        }

        let output_format = match self.output_format {
            Format::Auto => {
                if let Some(output_path) = self.output.as_deref() {
                    if output_path.ends_with(".sam") {
                        bam::Format::Sam
                    } else if output_path.ends_with(".cram") {
                        bam::Format::Cram
                    } else {
                        bam::Format::Bam
                    }
                } else {
                    bam::Format::Sam
                }
            }
            Format::BAM => bam::Format::Bam,
            Format::SAM => bam::Format::Sam,
            Format::CRAM => bam::Format::Cram,
        };

        let header_elements: Vec<_> = input_bam
            .header()
            .as_bytes()
            .split(|x| *x == b'\n')
            .map(|x| {
                if x.starts_with(b"@HD") {
                    let elements: Vec<_> = x
                        .split(|y| *y == b'\t')
                        .map(|y| {
                            if y.starts_with(b"SO:") {
                                b"SO:unsorted"
                            } else {
                                y
                            }
                        })
                        .collect();
                    elements.join(&b"\t"[..])
                } else {
                    x.to_vec()
                }
            })
            .collect();
        let mut header = bam::Header::new();
        for one in header_elements.iter() {
            if !one.is_empty() {
                header.push_record(&bam::header::HeaderRecord::new(&one[1..]));
            }
        }

        let mut output_bam = if let Some(bam_path) = self.output.as_deref() {
            bam::Writer::from_path(bam_path, &header, output_format)
        } else {
            bam::Writer::from_stdout(&header, output_format)
        }?;

        if let Some(reference) = self.reference.as_deref() {
            output_bam.set_reference(reference)?;
        }

        if let Some(write_threads) = self.write_threads {
            output_bam.set_threads(write_threads)?;
        }

        let mut pair_map: HashMap<Vec<u8>, bam::Record> = HashMap::new();
        let mut record = bam::Record::new();

        while let Some(result) = input_bam.read(&mut record) {
            match result {
                Ok(_) => {
                    if record.is_secondary() || record.is_supplementary() {
                        continue;
                    }
                    if let Some(paired_record) = pair_map.remove(record.qname()) {
                        if record.is_first_in_template() {
                            output_bam.write(&record)?;
                            output_bam.write(&paired_record)?;
                        } else {
                            output_bam.write(&paired_record)?;
                            output_bam.write(&record)?;
                        }
                    } else {
                        pair_map.insert(record.qname().to_vec(), record.clone());
                    }
                }
                Err(e) => return Err(e.into()),
            }
        }

        Ok(())
    }
}
