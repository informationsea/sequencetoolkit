use anyhow::Context;
use bio::io::fasta::IndexedReader;
use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, Record};
use rust_htslib::bam::HeaderView;
use rust_htslib::bcf::Read as _;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::io::prelude::*;
use std::str;

fn reverse_complement(acid: u8) -> u8 {
    match acid {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn reverse_complement_triplet(acid: [u8; 3]) -> [u8; 3] {
    [
        reverse_complement(acid[2]),
        reverse_complement(acid[1]),
        reverse_complement(acid[0]),
    ]
}

fn reverse_complement_vec(acid: &[u8]) -> Vec<u8> {
    let mut seq: Vec<_> = acid.iter().map(|x| reverse_complement(*x)).collect();
    seq.reverse();
    seq
}

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Mismatch {
    pub reference: u8,
    pub sequenced: u8,
}

impl Mismatch {
    pub fn new(reference: u8, sequenced: u8) -> Self {
        Mismatch {
            reference: reference.to_ascii_uppercase(),
            sequenced: sequenced.to_ascii_uppercase(),
        }
    }

    pub fn reverse_complement(&self) -> Mismatch {
        Mismatch {
            reference: reverse_complement(self.reference),
            sequenced: reverse_complement(self.sequenced),
        }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct MismatchTriplet {
    pub reference: [u8; 3],
    pub sequenced: [u8; 3],
}

impl MismatchTriplet {
    pub fn new(reference: [u8; 3], sequenced: [u8; 3]) -> Self {
        MismatchTriplet {
            reference: [
                reference[0].to_ascii_uppercase(),
                reference[1].to_ascii_uppercase(),
                reference[2].to_ascii_uppercase(),
            ],
            sequenced: [
                sequenced[0].to_ascii_uppercase(),
                sequenced[1].to_ascii_uppercase(),
                sequenced[2].to_ascii_uppercase(),
            ],
        }
    }

    pub fn reverse_complement(&self) -> MismatchTriplet {
        MismatchTriplet {
            reference: [
                reverse_complement(self.reference[2]),
                reverse_complement(self.reference[1]),
                reverse_complement(self.reference[0]),
            ],
            sequenced: [
                reverse_complement(self.sequenced[2]),
                reverse_complement(self.sequenced[1]),
                reverse_complement(self.sequenced[0]),
            ],
        }
    }
}

fn increment<T: PartialEq + Eq + Hash>(hash: &mut HashMap<T, usize>, key: T) {
    if let Some(v) = hash.get_mut(&key) {
        *v += 1;
    } else {
        hash.insert(key, 1);
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SequencingErrorCounts {
    pub total_reference_base: HashMap<u8, usize>,
    pub total_reference_triplet: HashMap<[u8; 3], usize>,
    pub total_sequenced_len: usize,
    pub total_reference_len: usize,
    pub mismatch: HashMap<Mismatch, usize>,
    pub mismatch_triplet: HashMap<MismatchTriplet, usize>,
    pub insertion_length: HashMap<u32, usize>,
    pub short_insertion: HashMap<Vec<u8>, usize>,
    pub deletion_length: HashMap<u32, usize>,
    pub short_deletion: HashMap<Vec<u8>, usize>,
    pub softclip_length: HashMap<u32, usize>,
    pub hardclip_length: HashMap<u32, usize>,
}

#[derive(Debug)]
pub struct SequencingErrorProcessor<R: Read + Seek> {
    count: SequencingErrorCounts,
    last_tid: u32,
    reference_seq: IndexedReader<R>,
    cached_reference_seq: Vec<u8>,
    cache_data: Vec<u8>,
    max_indel_length: usize,
    known_variants: Option<rust_htslib::bcf::IndexedReader>,
    known_variant_positions: HashSet<usize>,
    target_regions: Regions,
    name_to_len: HashMap<String, usize>,
}

impl<R: Read + Seek> SequencingErrorProcessor<R> {
    pub fn new(
        max_indel_length: usize,
        known_variants: Option<rust_htslib::bcf::IndexedReader>,
        reference_seq: IndexedReader<R>,
        regions: Regions,
    ) -> Self {
        let name_to_len: HashMap<_, _> = reference_seq
            .index
            .sequences()
            .iter()
            .map(|x| (x.name.to_string(), x.len as usize))
            .collect();

        SequencingErrorProcessor {
            count: SequencingErrorCounts::default(),
            last_tid: u32::MAX,
            reference_seq,
            cached_reference_seq: Vec::new(),
            cache_data: Vec::new(),
            max_indel_length,
            known_variants,
            known_variant_positions: HashSet::new(),
            target_regions: regions,
            name_to_len,
        }
    }

    pub fn count(&self) -> &SequencingErrorCounts {
        &self.count
    }

    pub fn load_known_variants(&mut self, chromosome: &[u8]) -> anyhow::Result<()> {
        self.known_variant_positions.clear();
        if let Some(reader) = self.known_variants.as_mut() {
            let rid = match reader.header().name2rid(chromosome) {
                Ok(rid) => rid,
                Err(rust_htslib::errors::Error::BcfUnknownContig { contig: _contig }) => {
                    return Ok(())
                }
                Err(e) => return Err(e).context("Failed to load rid of known variants"),
            };

            reader
                .fetch(rid, 0, None)
                .context("Failed to fetch known variants")?;
            let mut record = reader.empty_record();
            while let Some(r) = reader.read(&mut record) {
                r?;
                let ref_allele = record.alleles()[0];
                for one_alt in record.alleles().iter().skip(1) {
                    if ref_allele.len() > one_alt.len()
                        && ref_allele[0] == one_alt[0]
                        && one_alt.len() == 1
                    {
                        // DEL
                        // eprintln!(
                        //     "register deletion pos: {}",
                        //     TryInto::<usize>::try_into(record.pos()).unwrap() + 1
                        // );
                        for pos in (record.pos() + 1)..record.end() {
                            if self
                                .target_regions
                                .contains(chromosome, pos.try_into().unwrap())
                            {
                                self.known_variant_positions.insert(pos.try_into().unwrap());
                            }
                        }
                    } else if ref_allele.len() < one_alt.len()
                        && ref_allele[0] == one_alt[0]
                        && ref_allele.len() == 1
                    {
                        // INS
                        // eprintln!(
                        //     "register insertion pos: {}",
                        //     TryInto::<usize>::try_into(record.pos()).unwrap() + 1
                        // );
                        if self.target_regions.contains(
                            chromosome,
                            TryInto::<u64>::try_into(record.pos()).unwrap() + 1,
                        ) {
                            self.known_variant_positions
                                .insert(TryInto::<usize>::try_into(record.pos()).unwrap() + 1);
                        }
                    } else {
                        // SNV/MNV/Others
                        for pos in record.pos()..record.end() {
                            // eprintln!(
                            //     "register snv pos: {}",
                            //     TryInto::<usize>::try_into(record.pos()).unwrap() + 1
                            // );
                            if self
                                .target_regions
                                .contains(chromosome, pos.try_into().unwrap())
                            {
                                self.known_variant_positions.insert(pos.try_into().unwrap());
                            }
                        }
                    }
                }
            }
        }
        //eprintln!("known variants: {:?}", self.known_variant_positions);
        Ok(())
    }

    pub fn add_record(&mut self, header_view: &HeaderView, record: &Record) -> anyhow::Result<()> {
        self.cache_data.clear();
        let tid = if record.tid() >= 0 {
            record.tid().try_into().unwrap()
        } else {
            return Err(anyhow::anyhow!("Unmapped read: {}", record.tid()));
        };
        let seq_name = str::from_utf8(header_view.tid2name(tid))?;
        if !self
            .target_regions
            .regions
            .contains_key(seq_name.as_bytes())
        {
            return Ok(());
        }
        if tid != self.last_tid {
            self.reference_seq
                .fetch_all(seq_name)
                .with_context(|| format!("Failed to load {seq_name} from fastq"))?;
            self.cached_reference_seq.clear();
            self.reference_seq
                .read(&mut self.cached_reference_seq)
                .with_context(|| format!("Failed to read {seq_name} from fastq"))?;
            self.last_tid = tid;
            self.load_known_variants(seq_name.as_bytes())?;
        }

        let cigar_view = record.cigar();
        let reference_end: u64 = cigar_view.end_pos().try_into().unwrap();
        let seq_len: u64 = if let Some(l) = self.name_to_len.get(seq_name) {
            (*l).try_into().unwrap()
        } else {
            return Err(anyhow::anyhow!("Unknown sequence name: {}", seq_name));
        };

        let reversed = record.is_reverse();

        let reference_cache_start: usize = (0.max(record.pos() - 1)).try_into().unwrap();
        let reference_cache_end: usize = seq_len.min(reference_end + 1).try_into().unwrap();
        self.cache_data.extend_from_slice(
            &self.cached_reference_seq[reference_cache_start..reference_cache_end],
        );
        // reference.fetch(
        //     seq_name,
        //     (0.max(record.pos() - 1)).try_into().unwrap(),
        //     seq_len.min(reference_end + 1),
        // )?;
        // reference.read(&mut self.cache_data)?;
        if record.pos() == 0 {
            self.cache_data.insert(0, b'N');
        }
        if reference_end == seq_len {
            self.cache_data.push(b'N');
        }

        let seq = record.seq().as_bytes();

        self.count.total_sequenced_len += seq.len();
        self.count.total_reference_len += self.cache_data.len() - 2;
        //eprintln!("REF:{}", str::from_utf8(&self.cache_data)?);
        //eprintln!("ALT:{}", str::from_utf8(&seq)?);

        let mut seq_pos: usize = 0;
        let mut ref_pos: usize = 1;

        for one_cigar in cigar_view.iter() {
            match one_cigar {
                Cigar::Del(l) => {
                    let ref_pos_real = ref_pos + reference_cache_start;
                    // eprintln!("SEQUENCING ERROR: del ref pos real: {ref_pos_real}");
                    if !self.known_variant_positions.contains(&ref_pos_real)
                        && self
                            .target_regions
                            .contains(seq_name.as_bytes(), ref_pos_real.try_into().unwrap())
                    {
                        increment(&mut self.count.deletion_length, *l);
                        if (*l) as usize <= self.max_indel_length {
                            let mut del_seq: Vec<_> = self.cache_data
                                [ref_pos..(ref_pos + *l as usize)]
                                .iter()
                                .map(|x| x.to_ascii_uppercase())
                                .collect();
                            if reversed {
                                del_seq = reverse_complement_vec(&del_seq);
                            }
                            increment(&mut self.count.short_deletion, del_seq);
                        }
                    }
                    ref_pos += TryInto::<usize>::try_into(*l).unwrap();
                    //eprintln!("Del: {l}");
                }
                Cigar::Ins(l) => {
                    let ref_pos_real = ref_pos + reference_cache_start;
                    // eprintln!("SEQUENCING ERROR: ins ref pos real: {ref_pos_real}");
                    if !self.known_variant_positions.contains(&ref_pos_real)
                        && self
                            .target_regions
                            .contains(seq_name.as_bytes(), ref_pos_real.try_into().unwrap())
                    {
                        increment(&mut self.count.insertion_length, *l);
                        if (*l) as usize <= self.max_indel_length {
                            let mut ins_seq: Vec<_> = seq[seq_pos..(seq_pos + *l as usize)]
                                .to_vec()
                                .iter()
                                .map(|x| x.to_ascii_uppercase())
                                .collect();
                            if reversed {
                                ins_seq = reverse_complement_vec(&ins_seq);
                            }
                            increment(&mut self.count.short_insertion, ins_seq);
                        }
                    }
                    seq_pos += TryInto::<usize>::try_into(*l).unwrap();
                    //eprintln!("Ins: {l}");
                }
                Cigar::SoftClip(l) => {
                    let ref_pos_real = ref_pos + reference_cache_start;
                    if !self.known_variant_positions.contains(&ref_pos_real)
                        && self
                            .target_regions
                            .contains(seq_name.as_bytes(), ref_pos_real.try_into().unwrap())
                    {
                        increment(&mut self.count.softclip_length, *l);
                    }
                    seq_pos += TryInto::<usize>::try_into(*l).unwrap();
                    //eprintln!("Soft Clip: {l}");
                }
                Cigar::RefSkip(l) => {
                    ref_pos += TryInto::<usize>::try_into(*l).unwrap();
                }
                Cigar::HardClip(l) => {
                    let ref_pos_real = ref_pos + reference_cache_start;
                    if !self.known_variant_positions.contains(&ref_pos_real)
                        && self
                            .target_regions
                            .contains(seq_name.as_bytes(), ref_pos_real.try_into().unwrap())
                    {
                        increment(&mut self.count.hardclip_length, *l);
                    }
                }
                Cigar::Pad(_) => (),
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                    //eprintln!("Match: {l}");
                    let match_ref = &self.cache_data[ref_pos..(ref_pos + *l as usize)];
                    let match_seq = &seq[seq_pos..(seq_pos + *l as usize)];

                    //eprintln!("REF Match:{}", str::from_utf8(&match_ref)?);
                    //eprintln!("ALT Match:{}", str::from_utf8(&match_seq)?);

                    for (i, (r, s)) in match_ref.iter().zip(match_seq.iter()).enumerate() {
                        let ref_pos_real = ref_pos + reference_cache_start + i;
                        if !self.known_variant_positions.contains(&ref_pos_real)
                            && self
                                .target_regions
                                .contains(seq_name.as_bytes(), ref_pos_real.try_into().unwrap())
                        {
                            let mut ref_acid = (*r).to_ascii_uppercase();
                            if reversed {
                                ref_acid = reverse_complement(ref_acid)
                            }
                            increment(&mut self.count.total_reference_base, ref_acid);
                            let mut ref_triplet = [
                                self.cache_data[i + ref_pos - 1].to_ascii_uppercase(),
                                (*r).to_ascii_uppercase(),
                                self.cache_data[i + ref_pos + 1].to_ascii_uppercase(),
                            ];
                            if reversed {
                                ref_triplet = reverse_complement_triplet(ref_triplet);
                            }
                            increment(&mut self.count.total_reference_triplet, ref_triplet);
                            if (*r).to_ascii_uppercase() != (*s).to_ascii_uppercase() {
                                // eprintln!("SEQUENCING ERROR: snv/mnv ref pos real: {ref_pos_real}");
                                let mut mismatch = Mismatch::new(*r, *s);
                                let mut mismatch_triplet = MismatchTriplet::new(
                                    [
                                        self.cache_data[i + ref_pos - 1],
                                        *r,
                                        self.cache_data[i + ref_pos + 1],
                                    ],
                                    [
                                        self.cache_data[i + ref_pos - 1],
                                        *s,
                                        self.cache_data[i + ref_pos + 1],
                                    ],
                                );
                                if reversed {
                                    mismatch = mismatch.reverse_complement();
                                    mismatch_triplet = mismatch_triplet.reverse_complement();
                                }

                                increment(&mut self.count.mismatch, mismatch);
                                increment(&mut self.count.mismatch_triplet, mismatch_triplet);
                            }
                        }
                    }

                    seq_pos += TryInto::<usize>::try_into(*l).unwrap();
                    ref_pos += TryInto::<usize>::try_into(*l).unwrap();
                }
            }
        }
        Ok(())
    }

    pub fn add_bam(&mut self, bam_reader: &mut impl bam::Read, min_mapq: u8) -> anyhow::Result<()> {
        let header_view = bam_reader.header().clone();
        let mut record = Record::new();
        let mut record_count = 0;

        while let Some(r) = bam_reader.read(&mut record) {
            r.context("BAM parse error")?;
            if record.tid() < 0 {
                continue;
            }
            if record.mapq() < min_mapq {
                continue;
            }
            let seq_name = header_view.tid2name(record.tid().try_into().unwrap());
            if !self.target_regions.regions.contains_key(seq_name) {
                continue;
            }
            self.add_record(&header_view, &record)?;
            record_count += 1;
            if record_count % 10_000 == 0 {
                if record.seq_len() > 500 {
                    eprintln!("Processing {record_count} records");
                } else if record_count % 1_000_000 == 0 {
                    eprintln!("Processing {record_count} records");
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Regions {
    pub regions: HashMap<Vec<u8>, Vec<(u64, u64)>>,
}

impl Regions {
    pub fn new() -> Self {
        Regions {
            regions: HashMap::new(),
        }
    }

    pub fn create_from_fasta<R: Read + Seek>(
        reader: &mut bio::io::fasta::IndexedReader<R>,
    ) -> Self {
        let mut regions = Regions::new();
        for one_sequence in reader.index.sequences() {
            regions
                .regions
                .entry(one_sequence.name.as_bytes().to_vec())
                .or_insert_with(Vec::new)
                .push((0, one_sequence.len));
        }
        regions
    }

    pub fn load_from_bed(reader: impl Read) -> anyhow::Result<Self> {
        let mut regions = Regions::new();
        let mut bed_reader = bio::io::bed::Reader::new(reader);
        for r in bed_reader.records() {
            let r = r?;
            regions
                .regions
                .entry(r.chrom().as_bytes().to_vec())
                .or_insert_with(Vec::new)
                .push((r.start(), r.end()));
        }

        for v in regions.regions.values_mut() {
            v.sort();
            let mut merged_vec = Vec::new();
            let mut i = 0;
            while i < v.len() {
                let new_start = v[i].0;
                let mut new_end = v[i].1;
                let mut j = i + 1;
                while j < v.len() {
                    if v[j].0 <= new_end {
                        new_end = new_end.max(v[j].1);
                    } else {
                        break;
                    }
                    j += 1;
                }
                i = j;
                merged_vec.push((new_start, new_end));
            }
            v.clear();
            v.extend_from_slice(&merged_vec);
        }

        Ok(regions)
    }

    pub fn contains(&self, chromosome: &[u8], pos: u64) -> bool {
        if let Some(r) = self.regions.get(chromosome) {
            return r
                .binary_search_by(|probe| {
                    if probe.0 > pos {
                        std::cmp::Ordering::Greater
                    } else if probe.1 <= pos {
                        std::cmp::Ordering::Less
                    } else {
                        std::cmp::Ordering::Equal
                    }
                })
                .is_ok();
        }
        false
    }
}

#[cfg(test)]
mod test;
