use anyhow::Context;
use std::{
    io::Read,
    path::Path,
    process::{Child, Command},
};

#[derive(Debug)]
pub struct OraReader {
    ora_process: Child,
}

impl OraReader {
    pub fn new(
        orad_binary: impl AsRef<Path>,
        ora_fastq_path: impl AsRef<Path>,
    ) -> anyhow::Result<Self> {
        if !orad_binary.as_ref().exists() {
            return Err(anyhow::anyhow!(
                "Cannot find orad binary: {}",
                orad_binary.as_ref().display()
            ));
        }

        let ora_process = Command::new(orad_binary.as_ref())
            .arg("--stdout")
            .arg("--raw")
            .arg(ora_fastq_path.as_ref())
            .stdout(std::process::Stdio::piped())
            .spawn()
            .context("Cannot open FASTQ file with orad")?;
        Ok(OraReader { ora_process })
    }
}

impl Read for OraReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.ora_process.stdout.as_mut().unwrap().read(buf)
    }
}

pub fn open(
    orad_binary: impl AsRef<Path>,
    fastq_path: impl AsRef<Path>,
) -> anyhow::Result<Box<dyn Read + Send + Unpin>> {
    if let Some(ext) = fastq_path.as_ref().extension() {
        if ext == "ora" {
            return Ok(Box::new(OraReader::new(orad_binary, fastq_path)?));
        }
    }
    Ok(Box::new(autocompress::autodetect_open(fastq_path)?))
}
