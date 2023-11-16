use autocompress::io::{RayonReader, RayonWriter};
use clap::Args;
use std::io::{self, BufRead, Write};

#[derive(Debug, Args)]
#[command(
    about = "Remove non standard header and reorder header",
    version,
    author
)]
pub struct RemoveNonStandardHeader {
    #[arg(help = "Input VCF file")]
    input: Option<String>,
    #[arg(short, long, help = "Output VCF")]
    output: Option<String>,
    #[arg(short, long, help = "Add additional header line")]
    additional: Option<Vec<String>>,
}

impl RemoveNonStandardHeader {
    pub fn run(&self) -> anyhow::Result<()> {
        let mut reader =
            RayonReader::new(autocompress::autodetect_open_or_stdin(self.input.clone())?);
        let mut writer = RayonWriter::new(autocompress::autodetect_create_or_stdout_prefer_bgzip(
            self.input.as_deref(),
            autocompress::CompressionLevel::Default,
        )?);

        let mut header_file_format: Vec<Vec<u8>> = Vec::new();
        let mut header_info: Vec<Vec<u8>> = Vec::new();
        let mut header_format: Vec<Vec<u8>> = Vec::new();
        let mut header_alt: Vec<Vec<u8>> = Vec::new();
        let mut header_filter: Vec<Vec<u8>> = Vec::new();
        let mut header_contig: Vec<Vec<u8>> = Vec::new();

        let mut line: Vec<u8> = Vec::new();
        let mut index = 0;
        while reader.read_until(b'\n', &mut line)? > 0 {
            index += 1;
            if line.starts_with(b"##") {
                let header_line = vcf::VCFHeaderLine::from_bytes(&line, index)?;
                match header_line.contents() {
                    vcf::VCFHeaderContent::ALT { .. } => {
                        header_alt.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::INFO { .. } => {
                        header_info.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::FILTER { .. } => {
                        header_filter.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::FileFormat { .. } => {
                        header_file_format.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::Contig { .. } => {
                        header_contig.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::FORMAT { .. } => {
                        header_format.push(header_line.line().to_vec())
                    }
                    vcf::VCFHeaderContent::Other { .. } => (),
                }
            } else if line.starts_with(b"#") {
                if let Some(additional) = self.additional.as_ref() {
                    for one in header_file_format {
                        writer.write_all(&one)?;
                    }
                    for one in additional {
                        writer.write_all(one.as_bytes())?;
                        writer.write_all(b"\n")?;
                    }
                    for one_type in [
                        header_alt,
                        header_filter,
                        header_info,
                        header_format,
                        header_contig,
                    ]
                    .iter()
                    {
                        for one_line in one_type {
                            writer.write_all(one_line)?;
                        }
                    }
                }
                writer.write_all(&line)?;
                break;
            } else {
                writer.write_all(&line)?;
                break;
            }
            line.clear();
        }

        io::copy(&mut reader, &mut writer)?;

        Ok(())
    }
}
