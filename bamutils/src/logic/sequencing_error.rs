use bio::io::fasta::IndexedReader;
use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, Record};
use rust_htslib::bam::HeaderView;
use rust_htslib::bcf::Read as _;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::io::prelude::*;
use std::str;

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
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
}

fn increment<T: PartialEq + Eq + Hash>(hash: &mut HashMap<T, usize>, key: T) {
    if let Some(v) = hash.get_mut(&key) {
        *v += 1;
    } else {
        hash.insert(key, 1);
    }
}

#[derive(Default, Debug)]
pub struct SequencingErrorCount {
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

    last_tid: u32,
    cached_reference_seq: Vec<u8>,
    cache_data: Vec<u8>,
    max_indel_length: usize,
    known_variants: Option<rust_htslib::bcf::IndexedReader>,
    known_variant_positions: HashSet<usize>,
}

impl SequencingErrorCount {
    pub fn new(
        max_indel_length: usize,
        known_variants: Option<rust_htslib::bcf::IndexedReader>,
    ) -> Self {
        let mut s = SequencingErrorCount::default();
        s.max_indel_length = max_indel_length;
        s.last_tid = u32::MAX;
        s.known_variants = known_variants;
        s
    }

    pub fn load_known_variants(&mut self, chromosome: &[u8]) -> anyhow::Result<()> {
        self.known_variant_positions.clear();
        if let Some(reader) = self.known_variants.as_mut() {
            let rid = reader.header().name2rid(chromosome)?;
            reader.fetch(rid, 0, None)?;
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
                            self.known_variant_positions.insert(pos.try_into().unwrap());
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
                        self.known_variant_positions
                            .insert(TryInto::<usize>::try_into(record.pos()).unwrap() + 1);
                    } else {
                        // SNV/MNV/Others
                        for pos in record.pos()..record.end() {
                            // eprintln!(
                            //     "register snv pos: {}",
                            //     TryInto::<usize>::try_into(record.pos()).unwrap() + 1
                            // );
                            self.known_variant_positions.insert(pos.try_into().unwrap());
                        }
                    }
                }
            }
        }
        eprintln!("known variants: {:?}", self.known_variant_positions);
        Ok(())
    }

    pub fn add_record<R: Read + Seek>(
        &mut self,
        header_view: &HeaderView,
        record: &Record,
        reference: &mut IndexedReader<R>,
        name_to_len: &HashMap<String, usize>,
    ) -> anyhow::Result<()> {
        self.cache_data.clear();
        let tid = if record.tid() >= 0 {
            record.tid().try_into().unwrap()
        } else {
            return Err(anyhow::anyhow!("Unmapped read: {}", record.tid()));
        };
        let seq_name = str::from_utf8(header_view.tid2name(tid))?;
        if tid != self.last_tid {
            reference.fetch_all(seq_name)?;
            self.cached_reference_seq.clear();
            reference.read(&mut self.cached_reference_seq)?;
            self.last_tid = tid;
            self.load_known_variants(seq_name.as_bytes())?;
        }

        let cigar_view = record.cigar();
        let reference_end: u64 = cigar_view.end_pos().try_into().unwrap();
        let seq_len: u64 = if let Some(l) = name_to_len.get(seq_name) {
            (*l).try_into().unwrap()
        } else {
            return Err(anyhow::anyhow!("Unknown sequence name: {}", seq_name));
        };

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

        self.total_sequenced_len += seq.len();
        self.total_reference_len += self.cache_data.len() - 2;
        //eprintln!("REF:{}", str::from_utf8(&self.cache_data)?);
        //eprintln!("ALT:{}", str::from_utf8(&seq)?);

        let mut seq_pos: usize = 0;
        let mut ref_pos: usize = 1;

        for one_cigar in cigar_view.iter() {
            match one_cigar {
                Cigar::Del(l) => {
                    let ref_pos_real = ref_pos + reference_cache_start;
                    // eprintln!("SEQUENCING ERROR: del ref pos real: {ref_pos_real}");
                    if !self.known_variant_positions.contains(&ref_pos_real) {
                        increment(&mut self.deletion_length, *l);
                        if (*l) as usize <= self.max_indel_length {
                            increment(
                                &mut self.short_deletion,
                                self.cache_data[ref_pos..(ref_pos + *l as usize)]
                                    .iter()
                                    .map(|x| x.to_ascii_uppercase())
                                    .collect(),
                            );
                        }
                    }
                    ref_pos += TryInto::<usize>::try_into(*l).unwrap();
                    //eprintln!("Del: {l}");
                }
                Cigar::Ins(l) => {
                    let ref_pos_real = ref_pos + reference_cache_start;
                    // eprintln!("SEQUENCING ERROR: ins ref pos real: {ref_pos_real}");
                    if !self.known_variant_positions.contains(&ref_pos_real) {
                        increment(&mut self.insertion_length, *l);
                        if (*l) as usize <= self.max_indel_length {
                            increment(
                                &mut self.short_insertion,
                                seq[seq_pos..(seq_pos + *l as usize)]
                                    .to_vec()
                                    .iter()
                                    .map(|x| x.to_ascii_uppercase())
                                    .collect(),
                            );
                        }
                    }
                    seq_pos += TryInto::<usize>::try_into(*l).unwrap();
                    //eprintln!("Ins: {l}");
                }
                Cigar::SoftClip(l) => {
                    increment(&mut self.softclip_length, *l);
                    seq_pos += TryInto::<usize>::try_into(*l).unwrap();
                    //eprintln!("Soft Clip: {l}");
                }
                Cigar::RefSkip(l) => {
                    ref_pos += TryInto::<usize>::try_into(*l).unwrap();
                }
                Cigar::HardClip(l) => {
                    increment(&mut self.hardclip_length, *l);
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
                        if !self.known_variant_positions.contains(&ref_pos_real) {
                            increment(&mut self.total_reference_base, (*r).to_ascii_uppercase());
                            increment(
                                &mut self.total_reference_triplet,
                                [
                                    self.cache_data[i + ref_pos - 1].to_ascii_uppercase(),
                                    (*r).to_ascii_uppercase(),
                                    self.cache_data[i + ref_pos + 1].to_ascii_uppercase(),
                                ],
                            );
                            if (*r).to_ascii_uppercase() != (*s).to_ascii_uppercase() {
                                // eprintln!("SEQUENCING ERROR: snv/mnv ref pos real: {ref_pos_real}");
                                increment(&mut self.mismatch, Mismatch::new(*r, *s));
                                increment(
                                    &mut self.mismatch_triplet,
                                    MismatchTriplet::new(
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
                                    ),
                                );
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

    pub fn add_bam<R: Read + Seek>(
        &mut self,
        bam_reader: &mut impl bam::Read,
        reference: &mut IndexedReader<R>,
        min_mapq: u8,
    ) -> anyhow::Result<()> {
        let header_view = bam_reader.header().clone();
        let mut record = Record::new();
        let name_to_len: HashMap<_, _> = reference
            .index
            .sequences()
            .iter()
            .map(|x| (x.name.to_string(), x.len as usize))
            .collect();
        let mut record_count = 0;

        while let Some(r) = bam_reader.read(&mut record) {
            r?;
            if record.tid() < 0 {
                continue;
            }
            if record.mapq() < min_mapq {
                continue;
            }
            self.add_record(&header_view, &record, reference, &name_to_len)?;
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

#[cfg(test)]
mod test {
    use super::*;
    use rust_htslib::bam::Read as _;

    #[test]
    fn test_add_record() -> anyhow::Result<()> {
        let mut count = SequencingErrorCount::new(10, None);
        let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
        let name_to_len: HashMap<_, _> = reference_fasta
            .index
            .sequences()
            .iter()
            .map(|x| (x.name.to_string(), x.len as usize))
            .collect();
        let mut bam = bam::Reader::from_path("./testdata/demo1.bam")?;
        let mut record = bam::record::Record::new();

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ8");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert!(count.mismatch.is_empty());
        assert!(count.mismatch_triplet.is_empty());
        assert!(count.insertion_length.is_empty());
        assert!(count.deletion_length.is_empty());
        assert!(count.softclip_length.is_empty());
        assert_eq!(count.total_reference_len, 100);
        assert_eq!(count.total_sequenced_len, 100);

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ1");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert!(count.mismatch.is_empty());
        assert!(count.mismatch_triplet.is_empty());
        assert!(count.insertion_length.is_empty());
        assert!(count.deletion_length.is_empty());
        assert!(count.softclip_length.is_empty());
        assert_eq!(count.total_reference_len, 200);
        assert_eq!(count.total_sequenced_len, 200);

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ6");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert!(count.mismatch.is_empty());
        assert!(count.mismatch_triplet.is_empty());
        assert!(count.insertion_length.is_empty());
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 1), (2, 1), (3, 1)])
        );
        assert_eq!(
            count.short_deletion,
            HashMap::from([
                (b"A".to_vec(), 1),
                (b"AGC".to_vec(), 1),
                (b"CT".to_vec(), 1)
            ])
        );
        assert!(count.softclip_length.is_empty());

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ2");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 1),
                (Mismatch::new(b'T', b'A'), 1)
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TAT"), 1)
            ])
        );
        assert!(count.insertion_length.is_empty());
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 1), (2, 1), (3, 1)])
        );
        assert!(count.softclip_length.is_empty());

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ7");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 1),
                (Mismatch::new(b'T', b'A'), 1)
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TAT"), 1)
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 1), (2, 1), (3, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 1), (2, 1), (3, 1)])
        );

        assert_eq!(
            count.short_deletion,
            HashMap::from([
                (b"A".to_vec(), 1),
                (b"AGC".to_vec(), 1),
                (b"CT".to_vec(), 1)
            ])
        );

        assert_eq!(
            count.short_insertion,
            HashMap::from([
                (b"G".to_vec(), 1),
                (b"AA".to_vec(), 1),
                (b"TTT".to_vec(), 1)
            ])
        );
        assert!(count.softclip_length.is_empty());

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ5");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 2),
                (Mismatch::new(b'T', b'A'), 1),
                (Mismatch::new(b'G', b'T'), 1),
                (Mismatch::new(b'T', b'G'), 1)
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TAT"), 1),
                (MismatchTriplet::new(*b"GGC", *b"GTC"), 1),
                (MismatchTriplet::new(*b"CTC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TCC", *b"TGC"), 1),
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 1), (2, 1), (3, 1), (4, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 2), (2, 2), (3, 1), (4, 1)])
        );
        assert_eq!(count.softclip_length, HashMap::from([(9, 1)]));

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ3");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 2),
                (Mismatch::new(b'T', b'A'), 1),
                (Mismatch::new(b'G', b'T'), 1),
                (Mismatch::new(b'T', b'G'), 2),
                (Mismatch::new(b'A', b'G'), 2),
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TAT"), 1),
                (MismatchTriplet::new(*b"GGC", *b"GTC"), 1),
                (MismatchTriplet::new(*b"CTC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TCC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TAC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TGT"), 1),
                (MismatchTriplet::new(*b"AAA", *b"AGA"), 1),
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 1), (2, 1), (3, 1), (4, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 4), (2, 2), (3, 1), (4, 1)])
        );
        assert_eq!(count.softclip_length, HashMap::from([(9, 1)]));

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ4");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 2),
                (Mismatch::new(b'T', b'A'), 1),
                (Mismatch::new(b'G', b'T'), 1),
                (Mismatch::new(b'T', b'G'), 3),
                (Mismatch::new(b'A', b'G'), 2),
                (Mismatch::new(b'A', b'C'), 1),
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TAT"), 1),
                (MismatchTriplet::new(*b"GGC", *b"GTC"), 1),
                (MismatchTriplet::new(*b"CTC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TCC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TAC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TGT"), 1),
                (MismatchTriplet::new(*b"AAA", *b"AGA"), 1),
                (MismatchTriplet::new(*b"GAT", *b"GCT"), 1),
                (MismatchTriplet::new(*b"TTA", *b"TGA"), 1),
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 2), (2, 1), (3, 2), (4, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 4), (2, 2), (3, 1), (4, 1)])
        );
        assert_eq!(count.softclip_length, HashMap::from([(9, 1), (6, 1)]));

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ9");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;

        assert!(bam.read(&mut record).is_none());

        assert_eq!(count.total_sequenced_len, 900);
        assert_eq!(*count.total_reference_triplet.get(b"NGA").unwrap(), 1);
        assert_eq!(*count.total_reference_triplet.get(b"TGN").unwrap(), 1);

        Ok(())
    }

    #[test]
    fn test_add_bam() -> anyhow::Result<()> {
        let mut count = SequencingErrorCount::new(10, None);
        let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
        let mut bam = bam::Reader::from_path("./testdata/demo1.bam")?;
        count.add_bam(&mut bam, &mut reference_fasta, 0).unwrap();
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 2),
                (Mismatch::new(b'T', b'A'), 1),
                (Mismatch::new(b'G', b'T'), 1),
                (Mismatch::new(b'T', b'G'), 3),
                (Mismatch::new(b'A', b'G'), 2),
                (Mismatch::new(b'A', b'C'), 1),
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TAT"), 1),
                (MismatchTriplet::new(*b"GGC", *b"GTC"), 1),
                (MismatchTriplet::new(*b"CTC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TCC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TAC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TGT"), 1),
                (MismatchTriplet::new(*b"AAA", *b"AGA"), 1),
                (MismatchTriplet::new(*b"GAT", *b"GCT"), 1),
                (MismatchTriplet::new(*b"TTA", *b"TGA"), 1),
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 2), (2, 1), (3, 2), (4, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 4), (2, 2), (3, 1), (4, 1)])
        );
        assert_eq!(count.softclip_length, HashMap::from([(9, 1), (6, 1)]));
        Ok(())
    }

    #[test]
    fn test_add_cram() -> anyhow::Result<()> {
        let mut count = SequencingErrorCount::new(10, None);
        let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
        let mut bam = bam::Reader::from_path("./testdata/demo1.cram")?;
        bam.set_reference("./testdata/ref/MT.fa")?;
        count.add_bam(&mut bam, &mut reference_fasta, 0).unwrap();
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 2),
                (Mismatch::new(b'T', b'A'), 1),
                (Mismatch::new(b'G', b'T'), 1),
                (Mismatch::new(b'T', b'G'), 3),
                (Mismatch::new(b'A', b'G'), 2),
                (Mismatch::new(b'A', b'C'), 1),
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TAT"), 1),
                (MismatchTriplet::new(*b"GGC", *b"GTC"), 1),
                (MismatchTriplet::new(*b"CTC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TCC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TAC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TGT"), 1),
                (MismatchTriplet::new(*b"AAA", *b"AGA"), 1),
                (MismatchTriplet::new(*b"GAT", *b"GCT"), 1),
                (MismatchTriplet::new(*b"TTA", *b"TGA"), 1),
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 2), (2, 1), (3, 2), (4, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 4), (2, 2), (3, 1), (4, 1)])
        );
        assert_eq!(count.softclip_length, HashMap::from([(9, 1), (6, 1)]));
        Ok(())
    }

    #[test]
    fn test_add_record_with_known_variants() -> anyhow::Result<()> {
        let known_variants =
            rust_htslib::bcf::IndexedReader::from_path("./testdata/known_variants.vcf.gz")?;
        let mut count = SequencingErrorCount::new(10, Some(known_variants));
        let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;

        let name_to_len: HashMap<_, _> = reference_fasta
            .index
            .sequences()
            .iter()
            .map(|x| (x.name.to_string(), x.len as usize))
            .collect();
        let mut bam = bam::Reader::from_path("./testdata/demo1.bam")?;
        let mut record = bam::record::Record::new();

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ8");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert!(count.mismatch.is_empty());
        assert!(count.mismatch_triplet.is_empty());
        assert!(count.insertion_length.is_empty());
        assert!(count.deletion_length.is_empty());
        assert!(count.softclip_length.is_empty());
        assert_eq!(count.total_reference_len, 100);
        assert_eq!(count.total_sequenced_len, 100);

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ1");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert!(count.mismatch.is_empty());
        assert!(count.mismatch_triplet.is_empty());
        assert!(count.insertion_length.is_empty());
        assert!(count.deletion_length.is_empty());
        assert!(count.softclip_length.is_empty());
        assert_eq!(count.total_reference_len, 200);
        assert_eq!(count.total_sequenced_len, 200);

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ6");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert!(count.mismatch.is_empty());
        assert!(count.mismatch_triplet.is_empty());
        assert!(count.insertion_length.is_empty());
        assert_eq!(count.deletion_length, HashMap::from([(2, 1), (3, 1)]));
        assert_eq!(
            count.short_deletion,
            HashMap::from([(b"AGC".to_vec(), 1), (b"CT".to_vec(), 1)])
        );
        assert!(count.softclip_length.is_empty());

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ2");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([(Mismatch::new(b'C', b'G'), 1),])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([(MismatchTriplet::new(*b"CCC", *b"CGC"), 1),])
        );
        assert!(count.insertion_length.is_empty());
        assert_eq!(count.deletion_length, HashMap::from([(2, 1), (3, 1)]));
        assert!(count.softclip_length.is_empty());

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ7");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([(Mismatch::new(b'C', b'G'), 1),])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([(MismatchTriplet::new(*b"CCC", *b"CGC"), 1),])
        );
        assert_eq!(count.insertion_length, HashMap::from([(1, 1), (3, 1)]));
        assert_eq!(count.deletion_length, HashMap::from([(2, 1), (3, 1)]));

        assert_eq!(
            count.short_deletion,
            HashMap::from([(b"AGC".to_vec(), 1), (b"CT".to_vec(), 1)])
        );

        assert_eq!(
            count.short_insertion,
            HashMap::from([(b"G".to_vec(), 1), (b"TTT".to_vec(), 1)])
        );
        assert!(count.softclip_length.is_empty());

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ5");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 2),
                (Mismatch::new(b'G', b'T'), 1),
                (Mismatch::new(b'T', b'G'), 1)
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"GGC", *b"GTC"), 1),
                (MismatchTriplet::new(*b"CTC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TCC", *b"TGC"), 1),
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 1), (3, 1), (4, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 1), (2, 2), (3, 1), (4, 1)])
        );
        assert_eq!(count.softclip_length, HashMap::from([(9, 1)]));

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ3");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 2),
                (Mismatch::new(b'G', b'T'), 1),
                (Mismatch::new(b'T', b'G'), 2),
                (Mismatch::new(b'A', b'G'), 2),
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"GGC", *b"GTC"), 1),
                (MismatchTriplet::new(*b"CTC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TCC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TAC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TGT"), 1),
                (MismatchTriplet::new(*b"AAA", *b"AGA"), 1),
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 1), (3, 1), (4, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 3), (2, 2), (3, 1), (4, 1)])
        );
        assert_eq!(count.softclip_length, HashMap::from([(9, 1)]));

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ4");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;
        assert_eq!(
            count.mismatch,
            HashMap::from([
                (Mismatch::new(b'C', b'G'), 2),
                (Mismatch::new(b'G', b'T'), 1),
                (Mismatch::new(b'T', b'G'), 3),
                (Mismatch::new(b'A', b'G'), 2),
                (Mismatch::new(b'A', b'C'), 1),
            ])
        );
        assert_eq!(
            count.mismatch_triplet,
            HashMap::from([
                (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"GGC", *b"GTC"), 1),
                (MismatchTriplet::new(*b"CTC", *b"CGC"), 1),
                (MismatchTriplet::new(*b"TCC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TAC", *b"TGC"), 1),
                (MismatchTriplet::new(*b"TTT", *b"TGT"), 1),
                (MismatchTriplet::new(*b"AAA", *b"AGA"), 1),
                (MismatchTriplet::new(*b"GAT", *b"GCT"), 1),
                (MismatchTriplet::new(*b"TTA", *b"TGA"), 1),
            ])
        );
        assert_eq!(
            count.insertion_length,
            HashMap::from([(1, 2), (3, 2), (4, 1)])
        );
        assert_eq!(
            count.deletion_length,
            HashMap::from([(1, 3), (2, 2), (3, 1), (4, 1)])
        );
        assert_eq!(count.softclip_length, HashMap::from([(9, 1), (6, 1)]));

        bam.read(&mut record).unwrap()?;
        assert_eq!(record.qname(), b"READ9");
        count.add_record(bam.header(), &record, &mut reference_fasta, &name_to_len)?;

        assert!(bam.read(&mut record).is_none());

        assert_eq!(count.total_sequenced_len, 900);
        assert_eq!(*count.total_reference_triplet.get(b"NGA").unwrap(), 1);
        assert_eq!(*count.total_reference_triplet.get(b"TGN").unwrap(), 1);

        Ok(())
    }
}
