use super::*;
use crate::GeneAnnotError;
use log::warn;
use std::collections::{HashMap, HashSet};
use std::io::{self, BufReader};

pub fn load_refgene(
    genome: Genome,
    reader: impl io::Read,
) -> Result<GeneAnnotations, GeneAnnotError> {
    let table_reader = csv::ReaderBuilder::new()
        .escape(None)
        .quoting(false)
        .flexible(false)
        .delimiter(b'\t')
        .from_reader(BufReader::new(reader));
    let mut gene_to_transcript: HashMap<String, Vec<Transcript>> = HashMap::new();
    let mut warned_chromosomes = HashSet::new();

    for (line, record) in table_reader.into_records().enumerate() {
        let record = record?;
        if record.len() != 16 {
            return Err(GeneAnnotError::OtherError("missing columns").into());
        }
        let transcript_line = parse_transcript_line(&record)
            .map_err(|x| GeneAnnotError::RefGeneParseError(line as u64 + 1, Box::new(x)))?;
        if genome
            .chromosome_index(transcript_line.chromosome)
            .is_some()
        {
            let (gene_name, transcript) = parse_transcript(&transcript_line, &genome)
                .map_err(|x| GeneAnnotError::RefGeneParseError(line as u64 + 1, Box::new(x)))?;

            if let Some(l) = gene_to_transcript.get_mut(&gene_name) {
                l.push(transcript);
            } else {
                gene_to_transcript.insert(gene_name, vec![transcript]);
            }
        } else if !warned_chromosomes.contains(transcript_line.chromosome) {
            warn!("unknown chromosome name: {}", transcript_line.chromosome);
            warned_chromosomes.insert(transcript_line.chromosome.to_string());
        }
    }

    let genes: Vec<_> = gene_to_transcript
        .into_iter()
        .map(|(gene_name, transcripts)| Gene {
            id: gene_name.to_string(),
            symbol: gene_name,
            transcripts,
            annotations: HashMap::new(),
        })
        .collect();

    Ok(GeneAnnotations::new(genome, genes))
}

#[derive(Debug, PartialEq, Clone)]
struct TranscriptLine<'a> {
    transcript_id: &'a str,
    chromosome: &'a str,
    strand: Strand,
    tx_start: u64,
    tx_end: u64,
    cds_start: u64,
    cds_end: u64,
    exon_count: u64,
    exon_starts: Vec<u64>,
    exon_ends: Vec<u64>,
    score: u64,
    gene_id: &'a str,
    cds_start_stat: CdsStatus,
    cds_end_stat: CdsStatus,
    exon_frames: Vec<i8>,
}

fn parse_transcript_line(record: &csv::StringRecord) -> Result<TranscriptLine<'_>, GeneAnnotError> {
    Ok(TranscriptLine {
        transcript_id: record.get(1).unwrap(),
        chromosome: record.get(2).unwrap(),
        strand: match record.get(3).unwrap() {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unknown,
        },
        tx_start: record.get(4).unwrap().parse::<u64>()?,
        tx_end: record.get(5).unwrap().parse::<u64>()?,
        cds_start: record.get(6).unwrap().parse::<u64>()?,
        cds_end: record.get(7).unwrap().parse::<u64>()?,
        exon_count: record.get(8).unwrap().parse::<u64>()?,
        exon_starts: record
            .get(9)
            .unwrap()
            .split(',')
            .filter(|x| !x.is_empty())
            .try_fold::<_, _, Result<_, GeneAnnotError>>(Vec::new(), |mut v, x| {
                v.push(x.parse::<u64>()?);
                Ok(v)
            })?,
        exon_ends: record
            .get(10)
            .unwrap()
            .split(',')
            .filter(|x| !x.is_empty())
            .try_fold::<_, _, Result<_, GeneAnnotError>>(Vec::new(), |mut v, x| {
                v.push(x.parse::<u64>()?);
                Ok(v)
            })?,
        score: record.get(11).unwrap().parse::<u64>()?,
        gene_id: record.get(12).unwrap(),
        cds_start_stat: record.get(13).unwrap().parse()?,
        cds_end_stat: record.get(14).unwrap().parse()?,
        exon_frames: record
            .get(15)
            .unwrap()
            .split(',')
            .filter(|x| !x.is_empty())
            .try_fold::<_, _, Result<_, GeneAnnotError>>(Vec::new(), |mut v, x| {
                v.push(x.parse::<i8>()?);
                Ok(v)
            })?,
    })
}

fn parse_transcript(
    transcript_line: &TranscriptLine,
    genome: &Genome,
) -> Result<(String, Transcript), GeneAnnotError> {
    if transcript_line.cds_end == transcript_line.cds_start {
        let (name, transcript) = parse_noncoding_transcript(transcript_line, genome)?;
        Ok((name, Transcript::Noncoding(transcript)))
    } else {
        let (name, transcript) = parse_coding_transcript(transcript_line, genome)?;
        Ok((name, Transcript::Coding(transcript)))
    }
}

fn parse_coding_transcript(
    transcript_line: &TranscriptLine,
    genome: &Genome,
) -> Result<(String, CodingTranscript), GeneAnnotError> {
    let chromosome_index = genome.chromosome_index(transcript_line.chromosome).unwrap();
    let exons: Vec<_> = transcript_line
        .exon_starts
        .iter()
        .zip(transcript_line.exon_ends.iter())
        .map(|(s, e)| Exon {
            start: *s,
            end: *e,
            annotations: HashMap::new(),
        })
        .collect();

    Ok((
        transcript_line.gene_id.to_string(),
        CodingTranscript {
            id: transcript_line.transcript_id.to_string(),
            chromosome_index,
            strand: transcript_line.strand,
            start: transcript_line.tx_start,
            end: transcript_line.tx_end,
            cds_start: transcript_line.cds_start,
            cds_end: transcript_line.cds_end,
            cds_start_status: transcript_line.cds_start_stat,
            cds_end_status: transcript_line.cds_end_stat,
            annotations: HashMap::new(),
            exons,
        },
    ))
}

fn parse_noncoding_transcript(
    transcript_line: &TranscriptLine,
    genome: &Genome,
) -> Result<(String, NoncodingTranscript), GeneAnnotError> {
    let chromosome_index = genome.chromosome_index(transcript_line.chromosome).unwrap();
    Ok((
        transcript_line.gene_id.to_string(),
        NoncodingTranscript {
            id: transcript_line.transcript_id.to_string(),
            chromosome_index,
            strand: transcript_line.strand,
            start: transcript_line.tx_start,
            end: transcript_line.tx_end,
            annotations: HashMap::new(),
            exons: transcript_line
                .exon_starts
                .iter()
                .zip(transcript_line.exon_ends.iter())
                .map(|(s, e)| Exon {
                    start: *s,
                    end: *e,
                    annotations: HashMap::new(),
                })
                .collect(),
        },
    ))
}

#[cfg(test)]
mod test {
    use super::*;
    use bio::io::fasta::Index;

    #[allow(clippy::unreadable_literal)]
    #[test]
    fn test_parse_transcript_line() -> Result<(), GeneAnnotError> {
        let record = csv::StringRecord::from(vec![
            "777",
            "NM_001369787",
            "chr12",
            "-",
            "25205245",
            "25250929",
            "25209794",
            "25245384",
            "5",
            "25205245,25225613,25227233,25245273,25250763,",
            "25209911,25225773,25227412,25245395,25250929,",
            "0",
            "KRAS",
            "cmpl",
            "cmpl",
            "0,2,0,0,-1,",
        ]);
        assert_eq!(record.len(), 16);
        assert_eq!(
            parse_transcript_line(&record)?,
            TranscriptLine {
                transcript_id: "NM_001369787",
                chromosome: "chr12",
                strand: Strand::Reverse,
                tx_start: 25205245,
                tx_end: 25250929,
                cds_start: 25209794,
                cds_end: 25245384,
                exon_count: 5,
                exon_starts: vec![25205245, 25225613, 25227233, 25245273, 25250763],
                exon_ends: vec![25209911, 25225773, 25227412, 25245395, 25250929],
                score: 0,
                gene_id: "KRAS",
                cds_start_stat: CdsStatus::Complete,
                cds_end_stat: CdsStatus::Complete,
                exon_frames: vec![0, 2, 0, 0, -1],
            }
        );
        Ok(())
    }

    #[allow(clippy::unreadable_literal)]
    #[test]
    fn test_parse_coding_transcript() -> anyhow::Result<()> {
        let index =
            Index::from_file(&"testfiles/genome/human/GRCh38.primary_assembly.genome.fa.fai")?;
        let genome = fasta::load_fasta("GRCh38", &index);
        let record = csv::StringRecord::from(vec![
            "777",
            "NM_001369787",
            "chr12",
            "-",
            "25205245",
            "25250929",
            "25209794",
            "25245384",
            "5",
            "25205245,25225613,25227233,25245273,25250763,",
            "25209911,25225773,25227412,25245395,25250929,",
            "0",
            "KRAS",
            "cmpl",
            "cmpl",
            "0,2,0,0,-1,",
        ]);
        let transcript_line = parse_transcript_line(&record)?;

        assert_eq!(
            parse_coding_transcript(&transcript_line, &genome)?,
            (
                "KRAS".to_string(),
                CodingTranscript {
                    id: "NM_001369787".to_string(),
                    chromosome_index: 11,
                    strand: Strand::Reverse,
                    start: 25205245,
                    end: 25250929,
                    cds_start: 25209794,
                    cds_start_status: CdsStatus::Complete,
                    cds_end: 25245384,
                    cds_end_status: CdsStatus::Complete,
                    exons: vec![
                        Exon {
                            start: 25205245,
                            end: 25209911,
                            annotations: HashMap::new()
                        },
                        Exon {
                            start: 25225613,
                            end: 25225773,
                            annotations: HashMap::new()
                        },
                        Exon {
                            start: 25227233,
                            end: 25227412,
                            annotations: HashMap::new()
                        },
                        Exon {
                            start: 25245273,
                            end: 25245395,
                            annotations: HashMap::new()
                        },
                        Exon {
                            start: 25250763,
                            end: 25250929,
                            annotations: HashMap::new()
                        }
                    ],
                    annotations: HashMap::new()
                }
            )
        );

        Ok(())
    }

    #[allow(clippy::unreadable_literal)]
    #[test]
    fn test_parse_coding_transcript2() -> anyhow::Result<()> {
        let index =
            Index::from_file(&"testfiles/genome/human/GRCh38.primary_assembly.genome.fa.fai")?;
        let genome = fasta::load_fasta("GRCh38", &index);
        let record = csv::StringRecord::from(vec![
            "1559",
            "NM_001354870",
            "chr8",
            "+",
            "127735433",
            "127742951",
            "127736593",
            "127740958",
            "3",
            "127735433,127738250,127740395,",
            "127736623,127739019,127742951,",
            "0",
            "MYC",
            "cmpl",
            "cmpl",
            "0,0,1,",
        ]);
        let transcript_line = parse_transcript_line(&record)?;

        assert_eq!(
            parse_coding_transcript(&transcript_line, &genome)?,
            (
                "MYC".to_string(),
                CodingTranscript {
                    id: "NM_001354870".to_string(),
                    chromosome_index: 7,
                    strand: Strand::Forward,
                    start: 127735433,
                    end: 127742951,
                    cds_start: 127736593,
                    cds_start_status: CdsStatus::Complete,
                    cds_end: 127740958,
                    cds_end_status: CdsStatus::Complete,
                    exons: vec![
                        Exon {
                            start: 127735433,
                            end: 127736623,
                            annotations: HashMap::new()
                        },
                        Exon {
                            start: 127738250,
                            end: 127739019,
                            annotations: HashMap::new()
                        },
                        Exon {
                            start: 127740395,
                            end: 127742951,
                            annotations: HashMap::new()
                        },
                    ],
                    annotations: HashMap::new()
                }
            )
        );

        Ok(())
    }

    #[allow(clippy::unreadable_literal)]
    #[test]
    fn test_parse_noncoding_transcript() -> anyhow::Result<()> {
        let index =
            Index::from_file(&"testfiles/genome/human/GRCh38.primary_assembly.genome.fa.fai")?;
        let genome = fasta::load_fasta("GRCh38", &index);
        let record = csv::StringRecord::from(vec![
            "1444",
            "NR_164109",
            "chr13",
            "+",
            "112647043",
            "112684497",
            "112684497",
            "112684497",
            "3",
            "112647043,112679176,112683273,",
            "112647598,112679638,112684497,",
            "0",
            "ATP11AUN",
            "unk",
            "unk",
            "-1,-1,-1,",
        ]);
        let transcript_line = parse_transcript_line(&record)?;

        assert_eq!(
            parse_noncoding_transcript(&transcript_line, &genome)?,
            (
                "ATP11AUN".to_string(),
                NoncodingTranscript {
                    id: "NR_164109".to_string(),
                    chromosome_index: 12,
                    strand: Strand::Forward,
                    start: 112647043,
                    end: 112684497,
                    exons: vec![
                        Exon {
                            start: 112647043,
                            end: 112647598,
                            annotations: HashMap::new()
                        },
                        Exon {
                            start: 112679176,
                            end: 112679638,
                            annotations: HashMap::new()
                        },
                        Exon {
                            start: 112683273,
                            end: 112684497,
                            annotations: HashMap::new()
                        },
                    ],
                    annotations: HashMap::new()
                }
            )
        );

        Ok(())
    }

    use flate2::read::MultiGzDecoder;
    #[test]
    fn test_load_refgene() -> anyhow::Result<()> {
        let index = Index::from_file(&"testfiles/genome/human/hg38.fa.fai")?;
        let genome = fasta::load_fasta("GRCh38", &index);
        load_refgene(
            genome,
            MultiGzDecoder::new(&include_bytes!("../../../testfiles/gene/UCSC/refGene.txt.gz")[..]),
        )?;

        Ok(())
    }

    #[test]
    fn test_load_gencode() -> anyhow::Result<()> {
        let index = Index::from_file(&"testfiles/genome/human/hg38.fa.fai")?;
        let genome = fasta::load_fasta("GRCh38", &index);
        load_refgene(
            genome,
            MultiGzDecoder::new(
                &include_bytes!("../../../testfiles/gene/UCSC/wgEncodeGencodeBasicV33.txt.gz")[..],
            ),
        )?;

        Ok(())
    }
}
