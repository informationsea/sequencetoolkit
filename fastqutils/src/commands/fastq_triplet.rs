use anyhow::Result;
use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::process::{Command, Stdio};

#[derive(Parser, Debug, Clone, PartialEq)]
pub struct FastqTriplet {
    #[arg(help = "Input FASTQ (with gzip)")]
    fastq: String,
    #[arg(short, long, help = "Triplet output")]
    triplet_output: String,
    #[arg(short, long, help = "doublet output")]
    doublet_output: String,
    #[arg(short, long, help = "single output")]
    single_output: String,
}

impl FastqTriplet {
    pub fn run(&self) -> Result<()> {
        let mut gzip = Command::new("gzip");
        let mut p = gzip
            .arg("-dc")
            .arg(&self.fastq)
            .stdout(Stdio::piped())
            .spawn()
            .expect("Failed to start gzip");
        let mut stdout = BufReader::new(p.stdout.take().unwrap());
        let mut line = vec![];

        let mut single_map: HashMap<Vec<u8>, u64> = HashMap::new();
        let mut single_pos_map = HashMap::new();
        for i in 0..300 {
            single_pos_map.insert(i, HashMap::<Vec<u8>, u64>::new());
        }
        let mut doublet_map: HashMap<Vec<u8>, u64> = HashMap::new();
        let mut doublet_pos_map = HashMap::new();
        for i in 0..300 {
            doublet_pos_map.insert(i, HashMap::<Vec<u8>, u64>::new());
        }
        let mut triplet_map: HashMap<Vec<u8>, u64> = HashMap::new();
        let mut triplet_pos_map = HashMap::new();
        for i in 0..300 {
            triplet_pos_map.insert(i, HashMap::<Vec<u8>, u64>::new());
        }

        loop {
            // read header
            line.clear();
            let read_bytes = stdout.read_until(b'\n', &mut line)?;
            if read_bytes == 0 {
                break;
            }
            // read sequence
            line.clear();
            let read_bytes = stdout.read_until(b'\n', &mut line)?;
            if read_bytes == 0 {
                eprintln!("Unexpected EOF");
                break;
            }
            let sequence_line = line.strip_suffix(b"\n").unwrap();
            for i in 0..sequence_line.len() {
                let single = &sequence_line[i..i + 1];
                if let Some(v) = single_map.get_mut(single) {
                    *v += 1;
                } else {
                    single_map.insert(single.to_vec(), 1);
                }
                if let Some(v) = single_pos_map.get_mut(&i).unwrap().get_mut(single) {
                    *v += 1;
                } else {
                    single_pos_map
                        .get_mut(&i)
                        .unwrap()
                        .insert(single.to_vec(), 1);
                }
            }
            for i in 0..sequence_line.len() - 1 {
                let doublet = &sequence_line[i..i + 2];
                if let Some(v) = doublet_map.get_mut(doublet) {
                    *v += 1;
                } else {
                    doublet_map.insert(doublet.to_vec(), 1);
                }
                if let Some(v) = doublet_pos_map.get_mut(&i).unwrap().get_mut(doublet) {
                    *v += 1;
                } else {
                    doublet_pos_map
                        .get_mut(&i)
                        .unwrap()
                        .insert(doublet.to_vec(), 1);
                }
            }
            for i in 0..sequence_line.len() - 2 {
                let triplet = &sequence_line[i..i + 3];
                if let Some(v) = triplet_map.get_mut(triplet) {
                    *v += 1;
                } else {
                    triplet_map.insert(triplet.to_vec(), 1);
                }
                if let Some(v) = triplet_pos_map.get_mut(&i).unwrap().get_mut(triplet) {
                    *v += 1;
                } else {
                    triplet_pos_map
                        .get_mut(&i)
                        .unwrap()
                        .insert(triplet.to_vec(), 1);
                }
            }

            // read quality header
            line.clear();
            let read_bytes = stdout.read_until(b'\n', &mut line)?;
            if read_bytes == 0 {
                eprintln!("Unexpected EOF");
                break;
            }
            // read quality
            line.clear();
            let read_bytes = stdout.read_until(b'\n', &mut line)?;
            if read_bytes == 0 {
                eprintln!("Unexpected EOF");
                break;
            }
        }

        let mut triplet_output =
            File::create(&self.triplet_output).expect("Failed to open output file");
        for (k, v) in triplet_map.iter() {
            writeln!(triplet_output, "all\t{}\t{}", String::from_utf8_lossy(k), v)?;
        }
        for (k2, v2) in triplet_pos_map.iter() {
            for (k, v) in v2.iter() {
                writeln!(
                    triplet_output,
                    "{}\t{}\t{}",
                    k2,
                    String::from_utf8_lossy(k),
                    v
                )?;
            }
        }

        let mut doublet_output =
            File::create(&self.doublet_output).expect("Failed to open output file");
        for (k, v) in doublet_map.iter() {
            writeln!(doublet_output, "all\t{}\t{}", String::from_utf8_lossy(k), v)?;
        }
        for (k2, v2) in doublet_pos_map.iter() {
            for (k, v) in v2.iter() {
                writeln!(
                    doublet_output,
                    "{}\t{}\t{}",
                    k2,
                    String::from_utf8_lossy(k),
                    v
                )?;
            }
        }

        let mut single_output =
            File::create(&self.single_output).expect("Failed to open output file");
        for (k, v) in single_map.iter() {
            writeln!(single_output, "all\t{}\t{}", String::from_utf8_lossy(k), v)?;
        }
        for (k2, v2) in single_pos_map.iter() {
            for (k, v) in v2.iter() {
                writeln!(
                    single_output,
                    "{}\t{}\t{}",
                    k2,
                    String::from_utf8_lossy(k),
                    v
                )?;
            }
        }

        Ok(())
    }
}
