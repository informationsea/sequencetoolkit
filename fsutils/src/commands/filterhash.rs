use std::{
    fs::File,
    io::{Read, Write},
};

use anyhow::Context;
use clap::{Parser, ValueEnum};
use digest::{Digest, DynDigest};
use md5::Md5;
use sha1::Sha1;
use sha2::{Sha256, Sha512};

#[derive(ValueEnum, Debug, Clone, Copy, PartialEq)]
enum Algorithm {
    Sha1,
    Sha256,
    Sha512,
    Md5,
}

impl Algorithm {
    pub fn digest(&self) -> Box<dyn DynDigest> {
        match self {
            Algorithm::Sha1 => Box::new(Sha1::new()),
            Algorithm::Sha256 => Box::new(Sha256::new()),
            Algorithm::Sha512 => Box::new(Sha512::new()),
            Algorithm::Md5 => Box::new(Md5::new()),
        }
    }
}

#[derive(Parser, Debug, Clone, PartialEq)]
#[command(
    about = "Read input data from stdin and write output data to stdout while calculating hash",
    version,
    author
)]
pub struct FilterHash {
    #[arg(help = "Input files (read from stdin if not specified)")]
    input: Vec<String>,
    #[arg(short, long, help = "Output file (write to stdout if not specified)")]
    output: Option<String>,
    #[arg(
        short = 't',
        long,
        help = "Hash result output file (write to stderr if not specified)"
    )]
    hash_output: Option<String>,
    #[arg(help = "Hash Algorithm", short, long, default_value = "sha256")]
    algorithm: Algorithm,
}

impl FilterHash {
    pub fn run(&self) -> anyhow::Result<()> {
        let mut output_writer: Box<dyn Write> = self.output.as_ref().map_or_else(
            || -> std::io::Result<Box<dyn Write>> { Ok(Box::new(std::io::stdout())) },
            |x| File::create(x).map(|x| -> Box<dyn Write> { Box::new(x) }),
        )?;

        let mut hash_output_writer: Box<dyn Write> = self.hash_output.as_ref().map_or_else(
            || -> std::io::Result<Box<dyn Write>> { Ok(Box::new(std::io::stderr())) },
            |x| File::create(x).map(|x| -> Box<dyn Write> { Box::new(x) }),
        )?;

        let mut digest = self.algorithm.digest();

        let process_stdin = |output_writer: &mut Box<dyn Write>,
                             digest: &mut Box<dyn DynDigest>,
                             hash_output_writer: &mut Box<dyn Write>|
         -> anyhow::Result<()> {
            let hash = copy_and_hash(std::io::stdin().lock(), output_writer, digest)?;
            writeln!(hash_output_writer, "{}  -", hash)?;
            Ok(())
        };

        if self.input.is_empty() {
            process_stdin(&mut output_writer, &mut digest, &mut hash_output_writer)?;
        }

        for path in &self.input {
            if path == "-" {
                process_stdin(&mut output_writer, &mut digest, &mut hash_output_writer)?;
                continue;
            }
            let mut reader =
                File::open(path).with_context(|| format!("Failed to open {}", path))?;
            let hash = copy_and_hash(&mut reader, &mut output_writer, &mut digest)?;
            writeln!(hash_output_writer, "{}  {}", hash, path)?;
        }
        Ok(())
    }
}

fn copy_and_hash<R: Read, W: Write>(
    mut reader: R,
    mut writer: W,
    digest: &mut Box<dyn DynDigest>,
) -> anyhow::Result<String> {
    let mut buffer = [0; 1024 * 1024];
    digest.reset();
    loop {
        let read_size = reader.read(&mut buffer)?;
        if read_size == 0 {
            break;
        }
        digest.update(&buffer[..read_size]);
        writer.write_all(&buffer[..read_size])?;
    }
    Ok(digest
        .finalize_reset()
        .into_iter()
        .map(|x| format!("{:02x}", x))
        .collect())
}
