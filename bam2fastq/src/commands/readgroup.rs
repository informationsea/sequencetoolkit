use crate::utils::{aux_to_str, DNBSEQ_REGEX_IN_BAM, ILLUMINA_REGEX_IN_BAM};
use once_cell::sync::Lazy;
use regex::Regex;
use rust_htslib::bam::{self, Read};
use std::collections::HashSet;
use std::io::{BufWriter, Write};
use std::str;

#[derive(Debug, Clone, clap::Args)]
#[command(about = "print read groups", version, author)]
pub struct CountReadGroup {
    #[arg(help = "Input SAM/BAM/CRAM")]
    input: Option<String>,
    #[arg(long = "read-threads")]
    read_threads: Option<usize>,
    #[arg(long = "reference", short = 'T')]
    reference: Option<String>,
    #[arg(long = "do-not-print-header", short = 'h')]
    not_print_header: bool,
    #[arg(long = "output", short = 'o')]
    output: Option<String>,
}

pub static BWA_REGEX: Lazy<Regex> = Lazy::new(|| {
    regex::Regex::new(
        r"^@PG\tID:bwa.*CL:(?P<cmd>.*bwa mem.*-R\s(@RG[^\s]+ID:(?P<rg_id>[^\s\\]+)[^\s]+)(\s.+)?\s(?P<fastq1>[^\s]+)(\s+(?P<fastq2>[^\s]+))\s*)$",
    )
    .unwrap()
});

pub static BWA_MEM2_REGEX: Lazy<Regex> = Lazy::new(|| {
    regex::Regex::new(
        r"^@PG\tID:bwa-mem2.*CL:(?P<cmd>.*bwa-mem2 mem.*-R\s(@RG[^\s]+ID:(?P<rg_id>[^\s\\]+)[^\s]+)(\s.+)?\s(?P<fastq1>[^\s]+)(\s+(?P<fastq2>[^\s]+))\s*)$",
    )
    .unwrap()
});

impl CountReadGroup {
    pub fn run(&self) -> Result<(), anyhow::Error> {
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

        let mut output = BufWriter::new(autocompress::create_or_stdout(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?);

        if let Some(reference) = self.reference.as_deref() {
            input_bam.set_reference(reference)?;
        }

        if let Some(read_threads) = self.read_threads {
            input_bam.set_threads(read_threads)?;
        }
        let mut read_groups: Vec<_> = input_bam
            .header()
            .as_bytes()
            .split(|x| *x == b'\n')
            .filter(|x| x.starts_with(b"@RG"))
            .map(|x| ReadGroup {
                id: str::from_utf8(
                    &x.split(|y| *y == b'\t')
                        .filter(|y| y.starts_with(b"ID:"))
                        .next()
                        .expect("No read group ID")[3..],
                )
                .unwrap()
                .to_string(),
                sample_name: x
                    .split(|y| *y == b'\t')
                    .filter(|y| y.starts_with(b"SM:"))
                    .next()
                    .map(|y| str::from_utf8(&y[3..]).unwrap().to_string()),
                platform_unit: x
                    .split(|y| *y == b'\t')
                    .filter(|y| y.starts_with(b"PU:"))
                    .next()
                    .map(|y| str::from_utf8(&y[3..]).unwrap().to_string()),
                read_type: "".to_string(),
                read_prefix: "".to_string(),
                fastq1: "".to_string(),
                fastq2: "".to_string(),
                cmd: "".to_string(),
            })
            .collect();

        let mut unprocessed_read_group: HashSet<_> =
            read_groups.iter().map(|x| x.id.to_string()).collect();

        for record in input_bam.records() {
            if unprocessed_read_group.is_empty() {
                break;
            }

            let record = record?;
            let read_group = aux_to_str(&record.aux(b"RG")?).unwrap();
            if unprocessed_read_group.contains(read_group) {
                let read_name = str::from_utf8(record.qname()).unwrap();
                let (read_type, prefix) =
                    if let Some(cap) = ILLUMINA_REGEX_IN_BAM.captures(read_name) {
                        ("ILLUMINA", cap.name("prefix").unwrap().as_str())
                    } else if let Some(cap) = DNBSEQ_REGEX_IN_BAM.captures(read_name) {
                        ("DNBSEQ", cap.name("prefix").unwrap().as_str())
                    } else {
                        ("UNKNOWN", "")
                    };

                for one in read_groups.iter_mut() {
                    if one.id == read_group {
                        one.read_type = read_type.to_string();
                        one.read_prefix = prefix.to_string();
                    }
                }
                unprocessed_read_group.remove(read_group);
            }
        }

        let programs: Vec<_> = input_bam
            .header()
            .as_bytes()
            .split(|x| *x == b'\n')
            .filter(|x| x.starts_with(b"@PG"))
            .map(|x| str::from_utf8(x).unwrap())
            .filter_map(|x| BWA_REGEX.captures(x).or_else(|| BWA_MEM2_REGEX.captures(x)))
            .collect();
        for one_program in programs.iter() {
            let read_group = one_program.name("rg_id").unwrap().as_str();
            let fastq1 = one_program.name("fastq1").unwrap().as_str();
            let fastq2 = one_program.name("fastq2").unwrap().as_str();
            let cmd = one_program.name("cmd").unwrap().as_str();
            for one in read_groups.iter_mut() {
                if one.id == read_group {
                    one.fastq1 = fastq1.to_string();
                    one.fastq2 = fastq2.to_string();
                    one.cmd = cmd.to_string();
                }
            }
        }

        if !self.not_print_header {
            writeln!(output, "BAM_PATH\tREAD_GROUP\tSAMPLE_NAME\tPLATFORM_UNIT\tREAD_TYPE\tREAD_PREFIX\tFASTQ1\tFASTQ2\tCMD")?;
        }

        for one in read_groups {
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.input.as_deref().unwrap_or(""),
                one.id,
                one.sample_name.as_ref().map(|x| x.as_str()).unwrap_or(""),
                one.platform_unit.as_ref().map(|x| x.as_str()).unwrap_or(""),
                one.read_type,
                one.read_prefix,
                one.fastq1,
                one.fastq2,
                one.cmd
            )?;
        }
        Ok(())
    }
}

struct ReadGroup {
    id: String,
    sample_name: Option<String>,
    platform_unit: Option<String>,
    read_type: String,
    read_prefix: String,
    fastq1: String,
    fastq2: String,
    cmd: String,
}
