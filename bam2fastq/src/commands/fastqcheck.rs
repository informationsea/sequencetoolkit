use crate::utils::{DigestReader, NameType, DNBSEQ_REGEX, ILLUMINA_REGEX};
use anyhow::Context;
use autocompress::io::{RayonReader, RayonWriter};
use std::collections::HashMap;
use std::io::{prelude::*, BufReader};
use std::str;

#[derive(Debug, Clone, clap::Args)]
#[command(
    about = "Check order of FASTQ and create result summary",
    version,
    author
)]
pub struct FastqCheck {
    #[arg(help = "Input FASTQ 1")]
    fastq1: String,
    #[arg(help = "Input FASTQ 2")]
    fastq2: String,
    #[arg(help = "FASTQ check result", short = 'o', long = "output")]
    output: Option<String>,
    #[arg(
        help = "FASTQ check result summary",
        short = 's',
        long = "summary-output"
    )]
    summary_output: Option<String>,
    #[arg(
        help = "# of I/O threads",
        short = 't',
        long = "thread",
        default_value = "1"
    )]
    thread: usize,
    #[arg(
        short = 'b',
        long,
        help = "Path to orad binary",
        default_value = "orad"
    )]
    orad_binary: String,
}

const INITIAL_LOAD_ENTRIES: u64 = 5000;

impl FastqCheck {
    pub fn run(&self) -> Result<(), anyhow::Error> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.thread)
            .build_global()
            .context("Failed to set # of threads")?;

        let mut output = RayonWriter::new(autocompress::autodetect_create_or_stdout(
            self.output.as_deref(),
            autocompress::CompressionLevel::Default,
        )?);

        let summary_output = if let Some(summary_output) = self.summary_output.as_deref() {
            Some(RayonWriter::new(autocompress::autodetect_create(
                summary_output,
                autocompress::CompressionLevel::Default,
            )?))
        } else {
            None
        };

        writeln!(output, "FASTQ1\t{}", self.fastq1)?;
        writeln!(output, "FASTQ2\t{}", self.fastq2)?;

        match check_order(&self.orad_binary, &self.fastq1, &self.fastq2, &mut output) {
            Ok(summary) => {
                if let Some(mut summary_output) = summary_output {
                    writeln!(
                        summary_output,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\tOK",
                        self.output.as_deref().unwrap_or("stdout"),
                        self.fastq1,
                        self.fastq2,
                        summary.fastq_type,
                        summary.read_prefix,
                        summary.fastq1_md5,
                        summary.fastq2_md5,
                    )?;
                }
                drop(output);
                Ok(())
            }
            Err(e) => {
                if let Some(mut summary_output) = summary_output {
                    writeln!(
                        summary_output,
                        "{}\t{}\t{}\t\t\t\t\tNG",
                        self.output.as_deref().unwrap_or("stdout"),
                        self.fastq1,
                        self.fastq2,
                    )?;
                }
                writeln!(output, "ERROR\t{}", e)?;
                writeln!(output, "VERIFY_RESULT\tNG")?;
                output.flush()?;
                drop(output);
                Err(e)
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
struct FastqVerifyResult {
    fastq1_md5: String,
    fastq2_md5: String,
    read_prefix: String,
    fastq_type: NameType,
}

const RAYON_READER_CAPACITY: usize = 1024 * 1024 * 100;

fn check_order(
    orad_binary: impl AsRef<std::path::Path>,
    input_fastq1_path: impl AsRef<std::path::Path>,
    input_fastq2_path: impl AsRef<std::path::Path>,
    output: &mut impl Write,
) -> anyhow::Result<FastqVerifyResult> {
    if input_fastq1_path.as_ref() == input_fastq2_path.as_ref() {
        return Err(anyhow::anyhow!("FASTQ1 and FASTQ2 are same file."));
    }

    let mut input_fastq1_md5_reader = RayonReader::with_capacity(
        DigestReader::new(
            crate::utils::ora::open(orad_binary.as_ref(), input_fastq1_path)?,
            md5::Md5::default(),
        ),
        RAYON_READER_CAPACITY,
    );

    let mut input_fastq2_md5_reader = RayonReader::with_capacity(
        DigestReader::new(
            crate::utils::ora::open(orad_binary, input_fastq2_path)?,
            md5::Md5::default(),
        ),
        RAYON_READER_CAPACITY,
    );

    let mut input_fastq1 = FastqReadnameReader::new(BufReader::new(&mut input_fastq1_md5_reader));
    let mut input_fastq2 = FastqReadnameReader::new(BufReader::new(&mut input_fastq2_md5_reader));

    let mut initial_fastq1_readnames = Vec::new();
    let mut initial_fastq2_readnames = Vec::new();

    log::info!("start loading data");

    for _ in 0..INITIAL_LOAD_ENTRIES {
        let mut fastq1_line = String::new();
        let mut fastq2_line = String::new();

        if input_fastq1.next(&mut fastq1_line)? == 0 {
            if input_fastq2.next(&mut fastq2_line)? == 0 {
                break;
            }
            return Err(anyhow::anyhow!("Too small number of reads in FASTQ1"));
        }
        if input_fastq2.next(&mut fastq2_line)? == 0 {
            return Err(anyhow::anyhow!("Too small number of reads in FASTQ1"));
        }
        initial_fastq1_readnames.push(fastq1_line);
        initial_fastq2_readnames.push(fastq2_line);
    }

    let name_type = if initial_fastq1_readnames.is_empty() {
        NameType::Empty
    } else {
        NameType::suggest(&initial_fastq1_readnames[0], &initial_fastq2_readnames[0])?
    };

    log::info!("FASTQ type: {:?}", name_type);
    writeln!(output, "FASTQ_TYPE\t{}", name_type)?;

    let read_prefix = match name_type {
        NameType::Illumina => check_illumina_order(
            input_fastq1,
            input_fastq2,
            output,
            &initial_fastq1_readnames,
            &initial_fastq2_readnames,
        )?,
        NameType::DNBSeq => check_dnbseq_order(
            input_fastq1,
            input_fastq2,
            output,
            &initial_fastq1_readnames,
            &initial_fastq2_readnames,
        )?,
        NameType::Empty => "".to_string(),
    };

    log::info!("Finish FASTQ order check");

    let mut buffer = Vec::new();
    input_fastq1_md5_reader.read_to_end(&mut buffer)?;
    if !buffer.is_empty() {
        return Err(anyhow::anyhow!("FASTQ 1 reading is not completed."));
    }

    input_fastq2_md5_reader.read_to_end(&mut buffer)?;
    if !buffer.is_empty() {
        return Err(anyhow::anyhow!("FASTQ 2 reading is not completed."));
    }

    log::info!("Reading FASTQ completed");

    let input_fastq1_md5 = format!(
        "{:x}",
        input_fastq1_md5_reader.into_inner().finalize_reset()
    );
    let input_fastq2_md5 = format!(
        "{:x}",
        input_fastq2_md5_reader.into_inner().finalize_reset()
    );

    writeln!(output, "FASTQ1_MD5\t{}", input_fastq1_md5)?;
    writeln!(output, "FASTQ2_MD5\t{}", input_fastq2_md5)?;
    writeln!(output, "VERIFY_RESULT\tOK")?;

    log::info!("FASTQ MD5 {} {}", input_fastq1_md5, input_fastq2_md5);
    log::info!("Verify result OK");

    output.flush()?;

    Ok(FastqVerifyResult {
        fastq1_md5: input_fastq1_md5,
        fastq2_md5: input_fastq2_md5,
        read_prefix,
        fastq_type: name_type,
    })
}

fn check_dnbseq_order(
    mut input_fastq1: FastqReadnameReader<impl BufRead>,
    mut input_fastq2: FastqReadnameReader<impl BufRead>,
    output: &mut impl Write,
    initial_fastq1_readnames: &[String],
    initial_fastq2_readnames: &[String],
) -> anyhow::Result<String> {
    let cap1 = parse_dnbseq_readname(&initial_fastq1_readnames[0])?;
    writeln!(output, "PREFIX\t{}", cap1.prefix)?;

    let mut last_read1 = DNBSeqReadName {
        prefix: cap1.prefix,
        group: "",
        index: "",
        read: 0,
    };

    for (read1, read2) in initial_fastq1_readnames
        .iter()
        .zip(initial_fastq2_readnames.iter())
    {
        last_read1 = check_dnbseq_readname(&last_read1, output, read1, read2)?;
    }

    let mut fastq1_readname = String::new();
    let mut fastq2_readname = String::new();

    let mut fastq1_readname2 = String::new();
    let mut fastq2_readname2 = String::new();

    loop {
        let len1 = input_fastq1.next(&mut fastq1_readname)?;
        let len2 = input_fastq2.next(&mut fastq2_readname)?;

        if len1 == 0 && len2 == 0 {
            break;
        }

        last_read1 =
            check_dnbseq_readname(&last_read1, output, &fastq1_readname, &fastq2_readname)?;

        let len1 = input_fastq1.next(&mut fastq1_readname2)?;
        let len2 = input_fastq2.next(&mut fastq2_readname2)?;

        if len1 == 0 && len2 == 0 {
            break;
        }

        last_read1 =
            check_dnbseq_readname(&last_read1, output, &fastq1_readname2, &fastq2_readname2)?;
    }
    Ok(cap1.prefix.to_string())
}

fn check_illumina_order(
    mut input_fastq1: FastqReadnameReader<impl BufRead>,
    mut input_fastq2: FastqReadnameReader<impl BufRead>,
    mut output: &mut impl Write,
    initial_fastq1_readnames: &[String],
    initial_fastq2_readnames: &[String],
) -> anyhow::Result<String> {
    let cap1 = parse_illumina_readname(&initial_fastq1_readnames[0])?;
    let prefix = cap1.prefix.to_string();
    let mut processed_reads = initial_fastq1_readnames.len();

    // check most common index
    let mut index_occurrence = HashMap::new();
    for one in initial_fastq1_readnames {
        let cap1 = parse_illumina_readname(&one)?;
        //eprintln!("{:?}", cap1);
        if let Some(value) = index_occurrence.get_mut(cap1.read_index) {
            *value += 1;
        } else {
            index_occurrence.insert(cap1.read_index.to_string(), 1);
        }
    }
    let mut index_occurrence_vec: Vec<_> = index_occurrence.iter().collect();
    index_occurrence_vec.sort_by_key(|x| x.1);
    let common_index = index_occurrence_vec.last().unwrap().0;

    //eprintln!("fastq1 {}", initial_fastq1_readnames[0]);
    //eprintln!("fastq2 {}", initial_fastq2_readnames[0]);
    //eprintln!("common index: {}", common_index);
    log::debug!("occurrence: {:?}", index_occurrence_vec);

    writeln!(output, "READ_PREFIX\t{}", prefix)?;
    writeln!(output, "COMMON_INDEX\t{}", common_index)?;

    //writeln!(output, "SPECIAL_INDEX\tTILE\tX_POS\tY_POS\tREAD_INDEX")?;

    let mut fastq1_readname = String::new();
    let mut fastq2_readname = String::new();

    let mut last_tile_pos = (0, 0, 0);
    for (read1, read2) in initial_fastq1_readnames
        .iter()
        .zip(initial_fastq2_readnames.iter())
    {
        last_tile_pos = check_illumina_readname(
            &prefix,
            &common_index,
            last_tile_pos,
            &mut output,
            read1,
            read2,
        )?;
    }

    loop {
        if processed_reads % 10_000_000 == 0 {
            log::info!("Processed {} reads ", processed_reads,);
        }
        processed_reads += 1;
        let len1 = input_fastq1.next(&mut fastq1_readname)?;
        let len2 = input_fastq2.next(&mut fastq2_readname)?;

        if len1 == 0 && len2 == 0 {
            break;
        }

        last_tile_pos = check_illumina_readname(
            &prefix,
            &common_index,
            last_tile_pos,
            &mut output,
            &fastq1_readname,
            &fastq2_readname,
        )?;
    }

    Ok(prefix)
}

struct FastqReadnameReader<R: BufRead> {
    reader: R,
    line_index: u64,
    buffer: Vec<u8>,
}

impl<R: BufRead> FastqReadnameReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReadnameReader {
            reader,
            line_index: 0,
            buffer: Vec::new(),
        }
    }

    pub fn next(&mut self, buffer: &mut String) -> anyhow::Result<usize> {
        buffer.clear();
        let read_size = self.reader.read_line(buffer)?;
        self.line_index += 1;
        if read_size == 0 {
            return Ok(0);
        }
        if !buffer.starts_with('@') {
            return Err(anyhow::anyhow!("Bad FASTQ read name: {}", buffer));
        }

        self.buffer.clear();
        let size = self.reader.read_until(b'\n', &mut self.buffer)?;
        self.line_index += 1;
        if size == 0 {
            return Err(anyhow::anyhow!("Incomplete FASTQ entry: {}", buffer));
        }

        self.buffer.clear();
        let size = self.reader.read_until(b'\n', &mut self.buffer)?;
        self.line_index += 1;
        if size == 0 {
            return Err(anyhow::anyhow!("Incomplete FASTQ entry: {}", buffer));
        }

        self.buffer.clear();
        let size = self.reader.read_until(b'\n', &mut self.buffer)?;
        self.line_index += 1;
        if size == 0 {
            return Err(anyhow::anyhow!("Incomplete FASTQ entry: {}", buffer));
        }

        Ok(read_size)
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
struct IlluminaReadName<'a> {
    prefix: &'a str,
    tile: u16,
    y_pos: u16,
    x_pos: u16,
    read: u8,
    read_index: &'a str,
}

fn parse_illumina_readname(readname: &str) -> anyhow::Result<IlluminaReadName> {
    if let Some(cap) = ILLUMINA_REGEX.captures(readname) {
        Ok(IlluminaReadName {
            prefix: cap.name("prefix").unwrap().as_str(),
            tile: cap.name("tile").unwrap().as_str().parse().unwrap(),
            x_pos: cap.name("x_pos").unwrap().as_str().parse().unwrap(),
            y_pos: cap.name("y_pos").unwrap().as_str().parse().unwrap(),
            read: cap.name("read").unwrap().as_str().parse().unwrap(),
            read_index: cap.name("read_index").unwrap().as_str(),
        })
    } else {
        Err(anyhow::anyhow!("Not illumina FASTQ name: {}", readname))
    }
}

fn check_illumina_readname(
    expected_prefix: &str,
    common_index: &str,
    last_tile_pos: (u16, u16, u16),
    output: &mut impl Write,
    read1: &str,
    read2: &str,
) -> anyhow::Result<(u16, u16, u16)> {
    let cap1 = parse_illumina_readname(read1)?;
    let cap2 = parse_illumina_readname(read2)?;

    if cap1.prefix != expected_prefix {
        //eprintln!("Invalid prefix for read 1: {}", read1);
        return Err(anyhow::anyhow!("Invalid prefix for read 1: {}", read1));
    }

    if cap2.prefix != expected_prefix {
        //eprintln!("Invalid prefix for read 2: {}", read2);
        return Err(anyhow::anyhow!("Invalid prefix for read 2: {}", read2));
    }

    if cap1.tile != cap2.tile
        || cap1.x_pos != cap2.x_pos
        || cap1.y_pos != cap2.y_pos
        || cap1.read_index != cap2.read_index
    {
        return Err(anyhow::anyhow!(
            "Read 2 is not corresponds to read 1: {} / {}",
            read1,
            read2
        ));
    }

    if cap1.read != 1 {
        return Err(anyhow::anyhow!("Read 1 is not read 1: {}", read1));
    }

    if cap2.read != 2 {
        return Err(anyhow::anyhow!("Read 2 is not read 2: {}", read2));
    }

    if common_index != cap1.read_index {
        writeln!(
            output,
            "SPECIAL_INDEX\t{}\t{}\t{}\t{}",
            cap1.tile, cap1.x_pos, cap1.y_pos, cap1.read_index
        )?;
    }

    let next_tile_pos = (cap1.tile, cap1.y_pos, cap1.x_pos);
    if last_tile_pos >= next_tile_pos {
        return Err(anyhow::anyhow!(
            "Unexpected read order: {:?} - {:?} / {} / {}",
            last_tile_pos,
            next_tile_pos,
            read1,
            read2
        ));
    }

    Ok(next_tile_pos)
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
struct DNBSeqReadName<'a> {
    prefix: &'a str,
    group: &'a str,
    index: &'a str,
    read: u8,
}

fn parse_dnbseq_readname(readname: &str) -> anyhow::Result<DNBSeqReadName> {
    if let Some(cap) = DNBSEQ_REGEX.captures(readname) {
        Ok(DNBSeqReadName {
            prefix: cap.name("prefix").unwrap().as_str(),
            group: cap.name("group").unwrap().as_str(),
            index: cap.name("index").unwrap().as_str(),
            read: cap.name("read").unwrap().as_str().parse().unwrap(),
        })
    } else {
        Err(anyhow::anyhow!("Not DNBSeq FASTQ name: {}", readname))
    }
}

fn check_dnbseq_readname<'a>(
    last_name: &DNBSeqReadName,
    output: &mut impl Write,
    read1: &'a str,
    read2: &str,
) -> anyhow::Result<DNBSeqReadName<'a>> {
    let cap1 = parse_dnbseq_readname(read1)?;
    let cap2 = parse_dnbseq_readname(read2)?;

    if cap1.prefix != cap2.prefix
        || cap1.group != cap2.group
        || cap1.index != cap2.index
        || cap1.read != 1
        || cap2.read != 2
    {
        return Err(anyhow::anyhow!(
            "Read 2 is not corresponds to read 1: {} / {}",
            read1,
            read2
        ));
    }

    if last_name.prefix == cap1.prefix && last_name.group == cap1.group {
        if last_name >= &cap1 {
            return Err(anyhow::anyhow!(
                "Unexpected read order: {:?} -> {} / {}",
                last_name,
                read1,
                read2
            ));
        }
    } else {
        if !last_name.prefix.is_empty() && last_name.prefix != cap1.prefix {
            return Err(anyhow::anyhow!(
                "Unexpected flowcell: {:?} -> {} / {}",
                last_name,
                read1,
                read2
            ));
        }
        // new order
        writeln!(output, "GROUP_ORDER\t{}", cap1.group)?;
    }

    Ok(cap1)
}
