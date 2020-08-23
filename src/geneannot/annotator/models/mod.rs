use crate::geneannot::{GeneAnnotError, GeneAnnotErrorKind};
use bio::data_structures::interval_tree::IntervalTree;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::str::FromStr;

mod fasta;
mod refgene;

pub use fasta::load_fasta;
pub use refgene::load_refgene;
pub use serde_json::Value;

pub type Annotations = HashMap<String, Value>;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneAnnotations {
    genome: Genome,
    genes: Vec<Gene>,
    gene_to_index: HashMap<String, usize>,
    transcript_to_index: HashMap<String, (usize, usize)>,
    interval_tree: Vec<IntervalTree<u64, (usize, usize)>>,
}

impl GeneAnnotations {
    pub fn new(genome: Genome, genes: Vec<Gene>) -> GeneAnnotations {
        let gene_to_index: HashMap<_, _> = genes
            .iter()
            .enumerate()
            .map(|(i, x)| (x.id.to_string(), i))
            .collect();
        let transcript_to_index: HashMap<_, _> = genes
            .iter()
            .enumerate()
            .flat_map(|(i, x)| {
                x.transcripts()
                    .iter()
                    .enumerate()
                    .map(move |(j, y)| (y.id().to_string(), (i, j)))
            })
            .collect();
        let mut interval_tree: Vec<_> = genome
            .chromosomes()
            .iter()
            .map(|_| IntervalTree::new())
            .collect();
        for (i, gene) in genes.iter().enumerate() {
            for (j, transcript) in gene.transcripts().iter().enumerate() {
                interval_tree[transcript.chromosome_index()]
                    .insert(transcript.start()..transcript.end(), (i, j));
            }
        }

        GeneAnnotations {
            genome,
            genes,
            gene_to_index,
            transcript_to_index,
            interval_tree,
        }
    }
    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    pub fn genes(&self) -> &[Gene] {
        &self.genes
    }

    pub fn gene(&self, gene_id: &str) -> Option<&Gene> {
        self.gene_to_index.get(gene_id).map(|g| &self.genes[*g])
    }

    pub fn transcript(&self, transcript_id: &str) -> Option<(&Gene, &Transcript)> {
        self.transcript_to_index
            .get(transcript_id)
            .map(|(g, t)| (&self.genes[*g], &self.genes[*g].transcripts()[*t]))
    }
    pub fn interval_tree(
        &self,
        chromosome_index: usize,
    ) -> Option<&IntervalTree<u64, (usize, usize)>> {
        self.interval_tree.get(chromosome_index)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Genome {
    name: String,
    chromosomes: Vec<Chromosome>,
    name_to_index: HashMap<String, usize>,
}

impl Genome {
    pub fn new(name: &str, chromosomes: &[Chromosome]) -> Genome {
        Genome {
            name: name.to_string(),
            chromosomes: chromosomes.to_vec(),
            name_to_index: chromosomes
                .iter()
                .enumerate()
                .map(|(i, x)| (x.name.to_string(), i))
                .collect(),
        }
    }
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn chromosomes(&self) -> &[Chromosome] {
        &self.chromosomes
    }

    pub fn chromosome_index(&self, name: &str) -> Option<usize> {
        if let Some(index) = self.name_to_index.get(name) {
            Some(*index)
        } else if name.starts_with("chr") {
            self.name_to_index.get(&name[3..]).copied()
        } else {
            self.name_to_index.get(&format!("chr{}", name)).copied()
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Chromosome {
    pub name: String,
    pub length: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unknown => write!(f, "?"),
        }
    }
}

impl FromStr for Strand {
    type Err = GeneAnnotError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Reverse,
        })
    }
}

#[allow(clippy::len_without_is_empty)]
pub trait TranscriptTrait {
    fn id(&self) -> &str;
    fn exons(&self) -> &[Exon];
    fn chromosome_index(&self) -> usize;
    fn strand(&self) -> Strand;
    fn start(&self) -> u64;
    fn end(&self) -> u64;
    fn annotations(&self) -> &Annotations;
    fn cds_start(&self) -> Option<u64>;
    fn cds_end(&self) -> Option<u64>;
    fn cds_start_status(&self) -> Option<CdsStatus>;
    fn cds_end_status(&self) -> Option<CdsStatus>;
    fn cds_offset(&self) -> Option<u64> {
        match self.strand() {
            Strand::Forward | Strand::Unknown => self
                .cds_start()
                .map(|x| self.transcript_position(x).base_position()),
            Strand::Reverse => self
                .cds_end()
                .map(|x| self.transcript_position(x - 1).base_position()),
        }
    }
    fn cds_len(&self) -> Option<u64> {
        if let Some(cds_start) = self.cds_start() {
            if let Some(cds_end) = self.cds_end() {
                let cds_start_pos = self.transcript_position(cds_start).base_position();
                let cds_end_pos = self.transcript_position(cds_end).base_position();
                if cds_start_pos < cds_end_pos {
                    return Some(cds_end_pos - cds_start_pos);
                } else {
                    return Some(cds_start_pos - cds_end_pos);
                }
            }
        }
        None
    }
    fn len(&self) -> u64 {
        self.exons().iter().map(|x| x.len()).sum()
    }

    fn genome_position(&self, transcript_position: TranscriptPosition) -> u64 {
        let base_position = transcript_position.base_position();
        match self.strand() {
            Strand::Forward | Strand::Unknown => {
                let mut transcript_len = 0;
                let mut genome_position = self.start();
                for one_exon in self.exons() {
                    let exon_len = one_exon.len();
                    if base_position < transcript_len + exon_len {
                        genome_position = one_exon.start() + base_position - transcript_len;
                        break;
                    } else {
                        transcript_len += exon_len;
                    }
                }
                match transcript_position {
                    TranscriptPosition::Exon(_) => genome_position,
                    TranscriptPosition::AfterExon(_, x) => genome_position + x,
                    TranscriptPosition::BeforeExon(_, x) => genome_position - x,
                }
            }
            Strand::Reverse => {
                let mut transcript_len = 0;
                let mut genome_position = self.end();
                for one_exon in self.exons().iter().rev() {
                    let exon_len = one_exon.len();
                    if base_position < transcript_len + exon_len {
                        genome_position = one_exon.end() - (base_position - transcript_len) - 1;
                        break;
                    } else {
                        transcript_len += exon_len;
                    }
                }
                match transcript_position {
                    TranscriptPosition::Exon(_) => genome_position,
                    TranscriptPosition::AfterExon(_, x) => genome_position - x,
                    TranscriptPosition::BeforeExon(_, x) => genome_position + x,
                }
            }
        }
    }

    fn genome_position_from_cds(&self, position: CdsPosition) -> Option<u64> {
        if let Some(cds_start) = self.cds_start() {
            if let Some(cds_end) = self.cds_end() {
                let cds_offset = self
                    .transcript_position(match self.strand() {
                        Strand::Forward | Strand::Unknown => cds_start,
                        Strand::Reverse => cds_end - 1,
                    })
                    .base_position();
                let cds_end_offset = self
                    .transcript_position(match self.strand() {
                        Strand::Forward | Strand::Unknown => cds_end,
                        Strand::Reverse => cds_start - 1,
                    })
                    .base_position();
                let transcript_position = match position {
                    CdsPosition::BeforeCds(x) => {
                        x.with_new_base_position(cds_offset - x.base_position() - 1)
                    }
                    CdsPosition::Cds(x) => x.with_new_base_position(x.base_position() + cds_offset),
                    CdsPosition::AfterCds(x) => {
                        x.with_new_base_position(x.base_position() + cds_end_offset)
                    }
                };

                Some(self.genome_position(transcript_position))
            } else {
                None
            }
        } else {
            None
        }
    }

    fn cds_position(&self, position: u64) -> Option<CdsPosition> {
        if let Some(cds_start) = self.cds_start() {
            if let Some(cds_end) = self.cds_end() {
                let cds_offset = self
                    .transcript_position(match self.strand() {
                        Strand::Forward | Strand::Unknown => cds_start,
                        Strand::Reverse => cds_end - 1,
                    })
                    .base_position();
                let cds_end_offset = self
                    .transcript_position(match self.strand() {
                        Strand::Forward | Strand::Unknown => cds_end,
                        Strand::Reverse => cds_start - 1,
                    })
                    .base_position();

                let transcript_position = self.transcript_position(position);
                let base_position = transcript_position.base_position();
                if base_position < cds_offset {
                    Some(CdsPosition::BeforeCds(
                        transcript_position.with_new_base_position(cds_offset - base_position - 1),
                    ))
                } else if base_position < cds_end_offset {
                    Some(CdsPosition::Cds(
                        transcript_position.with_new_base_position(base_position - cds_offset),
                    ))
                } else {
                    Some(CdsPosition::AfterCds(
                        transcript_position.with_new_base_position(base_position - cds_end_offset),
                    ))
                }
            } else {
                None
            }
        } else {
            None
        }
    }

    fn transcript_position(&self, position: u64) -> TranscriptPosition {
        let exon_positions: Vec<_> = self
            .exons()
            .iter()
            .map(|x| x.exon_position(position, self.strand()))
            .collect();
        if exon_positions.iter().any(|x| match x {
            ExonPosition::InExon(_) => true,
            _ => false,
        }) {
            let pos = exon_positions
                .iter()
                .zip(self.exons().iter())
                .map(|(x, e)| match x {
                    ExonPosition::AfterExon(_) => e.end() - e.start(),
                    ExonPosition::InExon(p) => *p,
                    ExonPosition::BeforeExon(_) => 0,
                })
                .sum();
            TranscriptPosition::Exon(pos)
        } else {
            let nearest_exon = exon_positions
                .iter()
                .enumerate()
                .min_by_key(|x| x.1.distance())
                .unwrap();
            let length_sum_before_exons: u64 = exon_positions
                .iter()
                .zip(self.exons().iter())
                .map(|(x, e)| match x {
                    ExonPosition::AfterExon(_) => e.end() - e.start(),
                    _ => 0,
                })
                .sum();
            match nearest_exon.1 {
                ExonPosition::BeforeExon(x) => {
                    TranscriptPosition::BeforeExon(length_sum_before_exons, *x)
                }
                ExonPosition::AfterExon(x) => {
                    TranscriptPosition::AfterExon(length_sum_before_exons - 1, *x)
                }
                _ => unreachable!(),
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Transcript {
    Coding(CodingTranscript),
    Noncoding(NoncodingTranscript),
}

impl TranscriptTrait for Transcript {
    fn id(&self) -> &str {
        match self {
            Transcript::Coding(t) => t.id(),
            Transcript::Noncoding(t) => t.id(),
        }
    }
    fn exons(&self) -> &[Exon] {
        match self {
            Transcript::Coding(t) => t.exons(),
            Transcript::Noncoding(t) => t.exons(),
        }
    }
    fn chromosome_index(&self) -> usize {
        match self {
            Transcript::Coding(t) => t.chromosome_index(),
            Transcript::Noncoding(t) => t.chromosome_index(),
        }
    }
    fn strand(&self) -> Strand {
        match self {
            Transcript::Coding(t) => t.strand(),
            Transcript::Noncoding(t) => t.strand(),
        }
    }
    fn start(&self) -> u64 {
        match self {
            Transcript::Coding(t) => t.start(),
            Transcript::Noncoding(t) => t.start(),
        }
    }
    fn end(&self) -> u64 {
        match self {
            Transcript::Coding(t) => t.end(),
            Transcript::Noncoding(t) => t.end(),
        }
    }
    fn annotations(&self) -> &Annotations {
        match self {
            Transcript::Coding(t) => t.annotations(),
            Transcript::Noncoding(t) => t.annotations(),
        }
    }
    fn cds_start(&self) -> Option<u64> {
        match self {
            Transcript::Coding(t) => t.cds_start(),
            Transcript::Noncoding(t) => t.cds_start(),
        }
    }
    fn cds_end(&self) -> Option<u64> {
        match self {
            Transcript::Coding(t) => t.cds_end(),
            Transcript::Noncoding(t) => t.cds_end(),
        }
    }
    fn cds_start_status(&self) -> Option<CdsStatus> {
        match self {
            Transcript::Coding(t) => t.cds_start_status(),
            Transcript::Noncoding(t) => t.cds_start_status(),
        }
    }
    fn cds_end_status(&self) -> Option<CdsStatus> {
        match self {
            Transcript::Coding(t) => t.cds_end_status(),
            Transcript::Noncoding(t) => t.cds_end_status(),
        }
    }
    fn cds_offset(&self) -> Option<u64> {
        match self {
            Transcript::Coding(t) => t.cds_offset(),
            Transcript::Noncoding(t) => t.cds_offset(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CdsPosition {
    BeforeCds(TranscriptPosition),
    Cds(TranscriptPosition),
    AfterCds(TranscriptPosition),
}

impl PartialOrd for CdsPosition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CdsPosition {
    fn cmp(&self, other: &Self) -> Ordering {
        match self {
            CdsPosition::AfterCds(x) => match other {
                CdsPosition::AfterCds(y) => x.cmp(y),
                CdsPosition::Cds(_) => Ordering::Greater,
                CdsPosition::BeforeCds(_) => Ordering::Greater,
            },
            CdsPosition::Cds(x) => match other {
                CdsPosition::AfterCds(_) => Ordering::Less,
                CdsPosition::Cds(y) => x.cmp(y),
                CdsPosition::BeforeCds(_) => Ordering::Greater,
            },
            CdsPosition::BeforeCds(x) => match other {
                CdsPosition::AfterCds(_) => Ordering::Less,
                CdsPosition::Cds(_) => Ordering::Less,
                CdsPosition::BeforeCds(y) => x.cmp(y).reverse(),
            },
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TranscriptPosition {
    BeforeExon(u64, u64),
    Exon(u64),
    AfterExon(u64, u64),
}

impl PartialOrd for TranscriptPosition {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TranscriptPosition {
    fn cmp(&self, other: &Self) -> Ordering {
        self.base_position()
            .cmp(&other.base_position())
            .then_with(|| match self {
                TranscriptPosition::BeforeExon(_, x) => match other {
                    TranscriptPosition::BeforeExon(_, y) => x.cmp(y).reverse(),
                    TranscriptPosition::Exon(_) => Ordering::Less,
                    TranscriptPosition::AfterExon(_, _) => Ordering::Less,
                },
                TranscriptPosition::Exon(_) => match other {
                    TranscriptPosition::BeforeExon(_, _) => Ordering::Greater,
                    TranscriptPosition::Exon(_) => Ordering::Equal,
                    TranscriptPosition::AfterExon(_, _) => Ordering::Less,
                },
                TranscriptPosition::AfterExon(_, x) => match other {
                    TranscriptPosition::BeforeExon(_, _) => Ordering::Greater,
                    TranscriptPosition::Exon(_) => Ordering::Greater,
                    TranscriptPosition::AfterExon(_, y) => x.cmp(y),
                },
            })
    }
}

impl TranscriptPosition {
    pub fn base_position(&self) -> u64 {
        match self {
            TranscriptPosition::BeforeExon(x, _) => *x,
            TranscriptPosition::Exon(x) => *x,
            TranscriptPosition::AfterExon(x, _) => *x,
        }
    }

    pub fn with_new_base_position(&self, new_position: u64) -> TranscriptPosition {
        match self {
            TranscriptPosition::Exon(_) => TranscriptPosition::Exon(new_position),
            TranscriptPosition::BeforeExon(_, x) => {
                TranscriptPosition::BeforeExon(new_position, *x)
            }
            TranscriptPosition::AfterExon(_, x) => TranscriptPosition::AfterExon(new_position, *x),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum ExonPosition {
    BeforeExon(u64),
    InExon(u64),
    AfterExon(u64),
}

impl ExonPosition {
    pub fn distance(&self) -> u64 {
        match self {
            ExonPosition::BeforeExon(x) => *x,
            ExonPosition::InExon(_) => 0,
            ExonPosition::AfterExon(x) => *x,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Gene {
    id: String,
    symbol: String,
    transcripts: Vec<Transcript>,
    annotations: Annotations,
}

impl Gene {
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn symbol(&self) -> &str {
        &self.symbol
    }
    pub fn transcripts(&self) -> &[Transcript] {
        &self.transcripts
    }
    pub fn annotations(&self) -> &Annotations {
        &self.annotations
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum CdsStatus {
    Complete,
    Incomplete,
    Unknown,
    None,
}

impl FromStr for CdsStatus {
    type Err = GeneAnnotError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "none" => CdsStatus::None,
            "unk" => CdsStatus::Unknown,
            "incmpl" => CdsStatus::Incomplete,
            "cmpl" => CdsStatus::Complete,
            _ => {
                return Err(GeneAnnotErrorKind::OtherError("unknown CDS status").into());
            }
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CodingTranscript {
    id: String,
    exons: Vec<Exon>,
    chromosome_index: usize,
    strand: Strand,
    start: u64,
    end: u64,
    cds_start: u64,
    cds_end: u64,
    cds_start_status: CdsStatus,
    cds_end_status: CdsStatus,
    annotations: Annotations,
}

impl TranscriptTrait for CodingTranscript {
    fn id(&self) -> &str {
        &self.id
    }
    fn exons(&self) -> &[Exon] {
        &self.exons
    }
    fn chromosome_index(&self) -> usize {
        self.chromosome_index
    }
    fn strand(&self) -> Strand {
        self.strand
    }
    fn start(&self) -> u64 {
        self.start
    }
    fn end(&self) -> u64 {
        self.end
    }
    fn annotations(&self) -> &Annotations {
        &self.annotations
    }
    fn cds_start(&self) -> Option<u64> {
        Some(self.cds_start)
    }
    fn cds_end(&self) -> Option<u64> {
        Some(self.cds_end)
    }
    fn cds_start_status(&self) -> Option<CdsStatus> {
        Some(self.cds_start_status)
    }
    fn cds_end_status(&self) -> Option<CdsStatus> {
        Some(self.cds_end_status)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Exon {
    start: u64,
    end: u64,
    annotations: HashMap<String, String>,
}

#[allow(clippy::len_without_is_empty)]
impl Exon {
    pub fn start(&self) -> u64 {
        self.start
    }
    pub fn end(&self) -> u64 {
        self.end
    }
    pub fn annotations(&self) -> &HashMap<String, String> {
        &self.annotations
    }
    pub fn exon_position(&self, position: u64, strand: Strand) -> ExonPosition {
        match strand {
            Strand::Forward | Strand::Unknown => {
                if position < self.start {
                    ExonPosition::BeforeExon(self.start - position)
                } else if position < self.end {
                    ExonPosition::InExon(position - self.start)
                } else {
                    ExonPosition::AfterExon(position - self.end + 1)
                }
            }
            Strand::Reverse => {
                if position < self.start {
                    ExonPosition::AfterExon(self.start - position)
                } else if position < self.end {
                    ExonPosition::InExon(self.end - position - 1)
                } else {
                    ExonPosition::BeforeExon(position - self.end + 1)
                }
            }
        }
    }
    pub fn len(&self) -> u64 {
        self.end - self.start
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NoncodingTranscript {
    id: String,
    exons: Vec<Exon>,
    chromosome_index: usize,
    strand: Strand,
    start: u64,
    end: u64,
    annotations: Annotations,
}

impl TranscriptTrait for NoncodingTranscript {
    fn id(&self) -> &str {
        &self.id
    }
    fn exons(&self) -> &[Exon] {
        &self.exons
    }
    fn chromosome_index(&self) -> usize {
        self.chromosome_index
    }
    fn strand(&self) -> Strand {
        self.strand
    }
    fn start(&self) -> u64 {
        self.start
    }
    fn end(&self) -> u64 {
        self.end
    }
    fn annotations(&self) -> &Annotations {
        &self.annotations
    }
    fn cds_start(&self) -> Option<u64> {
        None
    }
    fn cds_end(&self) -> Option<u64> {
        None
    }
    fn cds_start_status(&self) -> Option<CdsStatus> {
        None
    }
    fn cds_end_status(&self) -> Option<CdsStatus> {
        None
    }
    fn cds_offset(&self) -> Option<u64> {
        None
    }
}

#[cfg(test)]
mod test;
