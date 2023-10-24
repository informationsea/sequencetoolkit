use flate2::write::GzEncoder;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::fs;
use std::io::{prelude::*, BufWriter};
use std::str;

#[derive(Debug, Clone, clap::Args)]
#[command(name = "bam2fastq", version, about = "Convert BAM/CRAM to FASTQ")]
pub struct Fastq {
    #[arg(help = "Input SAM/BAM/CRAM")]
    input: Option<String>,
    #[arg(short = 'o', long = "output-prefix")]
    output_prefix: String,
    #[arg(long = "read-threads")]
    read_threads: Option<usize>,
    #[arg(long = "write-threads", default_value = "1")]
    write_threads: usize,
    #[arg(long = "reference", short = 'r')]
    reference: Option<String>,
    #[arg(long = "compression-level", default_value = "6")]
    compression_level: u32,
}

impl Fastq {
    pub fn run(&self) -> anyhow::Result<()> {
        let mut input_bam = if let Some(bam_path) = self.input.as_deref() {
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

        let read_groups: Vec<_> = input_bam
            .header()
            .as_bytes()
            .split(|x| *x == b'\n')
            .filter(|x| x.starts_with(b"@RG"))
            .map(|x| {
                x.split(|y| *y == b'\t')
                    .filter(|y| y.starts_with(b"ID:"))
                    .next()
                    .expect("No read group ID")[3..]
                    .to_vec()
            })
            .collect();

        let iothread = autocompress::iothread::IoThread::new(self.write_threads);
        let compression_level = flate2::Compression::new(self.compression_level);

        let mut outputs = HashMap::new();
        for one_group in read_groups.iter() {
            outputs.insert(
                one_group.to_vec(),
                (
                    BufWriter::new(iothread.add_writer(GzEncoder::new(
                        fs::File::create(format!(
                            "{}_{}_1.fastq.gz",
                            self.output_prefix,
                            str::from_utf8(one_group).unwrap()
                        ))?,
                        compression_level,
                    ))),
                    BufWriter::new(iothread.add_writer(GzEncoder::new(
                        fs::File::create(format!(
                            "{}_{}_2.fastq.gz",
                            self.output_prefix,
                            str::from_utf8(one_group).unwrap()
                        ))?,
                        compression_level,
                    ))),
                ),
            );
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
                        let read_group =
                            crate::utils::aux_to_string(&record.aux(b"RG")?).expect("Invalid RG");
                        let writers = outputs.get_mut(read_group).ok_or_else(|| {
                            anyhow::anyhow!("Unknown RG: {}", str::from_utf8(read_group).unwrap())
                        })?;
                        if record.is_first_in_template() {
                            write_record(&mut writers.0, &record)?;
                            write_record(&mut writers.1, &paired_record)?;
                        } else {
                            write_record(&mut writers.1, &record)?;
                            write_record(&mut writers.0, &paired_record)?;
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

fn write_record<W: Write>(mut writer: W, record: &bam::Record) -> anyhow::Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(record.qname())?;
    writer.write_all(b"\n")?;

    if record.is_reverse() {
        let mut seq: Vec<_> = record
            .seq()
            .as_bytes()
            .iter()
            .map(|x| match x {
                b'A' => b'T',
                b'C' => b'G',
                b'T' => b'A',
                b'G' => b'C',
                _ => b'N',
            })
            .collect();
        seq.reverse();
        writer.write_all(&seq)?;
    } else {
        writer.write_all(&record.seq().as_bytes())?;
    }

    writer.write_all(b"\n+\n")?;
    let mut qual: Vec<_> = record.qual().iter().map(|x| x + 33).collect();
    if record.is_reverse() {
        qual.reverse();
    }
    writer.write_all(&qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}
