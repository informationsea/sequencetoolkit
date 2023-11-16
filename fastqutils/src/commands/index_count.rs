use anyhow::Context;
use clap::Parser;
use std::io::prelude::*;
use std::{collections::HashMap, io::BufRead};

#[derive(Parser, Debug, Clone, PartialEq)]
#[command(about = "Count indexes", version, author)]
pub struct IndexCount {
    #[arg(help = "Input FASTQ Read 1 (with gzip)")]
    fastq: String,
    #[arg(short, long, help = "output")]
    output: String,
    #[arg(short, long, help = "# of threads", default_value = "1")]
    threads: usize,
}

impl IndexCount {
    pub fn run(&self) -> anyhow::Result<()> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .context("Failed to set # of threads")?;
        let fastq = std::io::BufReader::new(autocompress::autodetect_open(&self.fastq)?);
        let mut output = std::io::BufWriter::new(std::fs::File::create(&self.output)?);

        let count = count_index(fastq)?;

        let mut fastq1_counts: Vec<_> = count.fastq1.iter().collect();
        fastq1_counts.sort_by_key(|x| x.1);
        fastq1_counts.reverse();
        let mut fastq2_counts: Vec<_> = count.fastq2.iter().collect();
        fastq2_counts.sort_by_key(|x| x.1);
        fastq2_counts.reverse();
        let mut pair_counts: Vec<_> = count.pair.iter().collect();
        pair_counts.sort_by_key(|x| x.1);
        pair_counts.reverse();

        for (index, count) in fastq1_counts {
            writeln!(
                &mut output,
                "read1\t{}\t{}",
                String::from_utf8_lossy(index),
                count
            )?;
        }

        for (index, count) in fastq2_counts {
            writeln!(
                output,
                "read2\t{}\t{}",
                String::from_utf8_lossy(index),
                count
            )?;
        }

        for (index, count) in pair_counts {
            writeln!(
                output,
                "pair\t{}\t{}",
                String::from_utf8_lossy(&index),
                count
            )?;
        }

        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
struct Count {
    fastq1: HashMap<Vec<u8>, u64>,
    fastq2: HashMap<Vec<u8>, u64>,
    pair: HashMap<Vec<u8>, u64>,
}

fn count_index<R1: BufRead>(mut fastq: R1) -> anyhow::Result<Count> {
    let mut count = Count::default();

    let mut line_count = 0;
    let mut buffer = Vec::new();
    loop {
        buffer.clear();

        let buf_size = fastq.read_until(b'\n', &mut buffer)?;

        if buf_size == 0 {
            break;
        }

        match line_count % 4 {
            0 => {
                // read name
                let index = String::from_utf8_lossy(
                    buffer
                        .split(|x| *x == b':')
                        .last()
                        .expect("No index found 1"),
                )
                .trim()
                .to_string();
                let indexes: Vec<_> = index.split('+').map(|x| x.as_bytes().to_vec()).collect();
                if indexes.len() != 2 {
                    anyhow::bail!("Invalid index: {}", index);
                }

                if let Some(val) = count.pair.get_mut(index.as_bytes()) {
                    *val += 1;
                } else {
                    count.pair.insert(index.as_bytes().to_vec(), 1);
                }

                if let Some(val) = count.fastq1.get_mut(&indexes[0]) {
                    *val += 1;
                } else {
                    count.fastq1.insert(indexes[0].to_vec(), 1);
                }

                if let Some(val) = count.fastq2.get_mut(&indexes[1]) {
                    *val += 1;
                } else {
                    count.fastq2.insert(indexes[1].to_vec(), 1);
                }
            }
            1 => {
                // sequence
            }
            2 => {
                // read name
            }
            3 => {
                // quality
            }
            _ => {
                unreachable!("line_count = {}", line_count)
            }
        }

        line_count += 1;
    }

    Ok(count)
}
