use bio::io::fasta::IndexedReader;
use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, Record};
use rust_htslib::bam::HeaderView;
use std::collections::HashMap;
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
            reference,
            sequenced,
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
            reference,
            sequenced,
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

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct SequencingErrorCount {
    pub total_match_base: HashMap<u8, usize>,
    pub total_match_triplet: HashMap<[u8; 3], usize>,
    pub total_sequenced_len: usize,
    pub total_reference_len: usize,
    pub mismatch: HashMap<Mismatch, usize>,
    pub mismatch_triplet: HashMap<MismatchTriplet, usize>,
    pub insertion_length: HashMap<u32, usize>,
    pub deletion_length: HashMap<u32, usize>,
    pub softclip_length: HashMap<u32, usize>,
    pub hardclip_length: HashMap<u32, usize>,

    cache_data: Vec<u8>,
}

impl SequencingErrorCount {
    pub fn new() -> Self {
        SequencingErrorCount::default()
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
        let cigar_view = record.cigar();
        let reference_end: u64 = cigar_view.end_pos().try_into().unwrap();
        let seq_len: u64 = if let Some(l) = name_to_len.get(seq_name) {
            (*l).try_into().unwrap()
        } else {
            return Err(anyhow::anyhow!("Unknown sequence name: {}", seq_name));
        };

        reference.fetch(
            seq_name,
            (0.max(record.pos() - 1)).try_into().unwrap(),
            seq_len.min(reference_end + 1),
        )?;
        reference.read(&mut self.cache_data)?;
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
                    increment(&mut self.deletion_length, *l);
                    ref_pos += (*l) as usize;
                    //eprintln!("Del: {l}");
                }
                Cigar::Ins(l) => {
                    increment(&mut self.insertion_length, *l);
                    seq_pos += (*l) as usize;
                    //eprintln!("Ins: {l}");
                }
                Cigar::SoftClip(l) => {
                    increment(&mut self.softclip_length, *l);
                    seq_pos += (*l) as usize;
                    //eprintln!("Soft Clip: {l}");
                }
                Cigar::RefSkip(l) => {
                    ref_pos += (*l) as usize;
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
                        increment(&mut self.total_match_base, *r);
                        increment(
                            &mut self.total_match_triplet,
                            [
                                self.cache_data[i + ref_pos - 1],
                                *r,
                                self.cache_data[i + ref_pos + 1],
                            ],
                        );
                        if *r != *s {
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

                    seq_pos += (*l) as usize;
                    ref_pos += (*l) as usize;
                }
            }
        }
        Ok(())
    }

    pub fn add_bam<R: Read + Seek>(
        &mut self,
        bam_reader: &mut impl bam::Read,
        reference: &mut IndexedReader<R>,
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
            self.add_record(&header_view, &record, reference, &name_to_len)?;
            record_count += 1;
            if record_count % 10_000 == 0 {
                if record.seq_len() > 500 {
                    eprintln!("Processing {record_count} records");
                } else if record_count % 100_000 == 0 {
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
        let mut count = SequencingErrorCount::new();
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
        assert_eq!(*count.total_match_triplet.get(b"NGA").unwrap(), 1);
        assert_eq!(*count.total_match_triplet.get(b"TGN").unwrap(), 1);

        Ok(())
    }

    #[test]
    fn test_add_bam() -> anyhow::Result<()> {
        let mut count = SequencingErrorCount::new();
        let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
        let mut bam = bam::Reader::from_path("./testdata/demo1.bam")?;
        count.add_bam(&mut bam, &mut reference_fasta).unwrap();
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
        let mut count = SequencingErrorCount::new();
        let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
        let mut bam = bam::Reader::from_path("./testdata/demo1.cram")?;
        bam.set_reference("./testdata/ref/MT.fa")?;
        count.add_bam(&mut bam, &mut reference_fasta).unwrap();
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
}
