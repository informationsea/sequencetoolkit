use std::fs::File;

use super::*;
use rust_htslib::bam::Read as _;

#[test]
fn test_add_record() -> anyhow::Result<()> {
    let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
    let regions = Regions::create_from_fasta(&mut reference_fasta);
    let mut processor = SequencingErrorProcessor::new(10, None, reference_fasta, regions);

    let mut bam = bam::Reader::from_path("./testdata/demo1.bam")?;
    let mut record = bam::record::Record::new();

    bam.read(&mut record).unwrap()?;
    assert_eq!(record.qname(), b"READ8");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
    assert!(count.mismatch.is_empty());
    assert!(count.mismatch_triplet.is_empty());
    assert!(count.insertion_length.is_empty());
    assert!(count.deletion_length.is_empty());
    assert!(count.softclip_length.is_empty());
    assert_eq!(count.total_reference_len, 100);
    assert_eq!(count.total_sequenced_len, 100);

    bam.read(&mut record).unwrap()?;
    assert_eq!(record.qname(), b"READ1");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
    assert!(count.mismatch.is_empty());
    assert!(count.mismatch_triplet.is_empty());
    assert!(count.insertion_length.is_empty());
    assert!(count.deletion_length.is_empty());
    assert!(count.softclip_length.is_empty());
    assert_eq!(count.total_reference_len, 200);
    assert_eq!(count.total_sequenced_len, 200);

    bam.read(&mut record).unwrap()?;
    assert_eq!(record.qname(), b"READ6");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    assert_eq!(record.qname(), b"READ10");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
    assert_eq!(
        count.mismatch,
        HashMap::from([
            (Mismatch::new(b'C', b'G'), 2),
            (Mismatch::new(b'T', b'A'), 2),
            (Mismatch::new(b'G', b'T'), 1),
            (Mismatch::new(b'T', b'G'), 3),
            (Mismatch::new(b'A', b'G'), 2),
            (Mismatch::new(b'A', b'C'), 3),
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
            (MismatchTriplet::new(*b"TTG", *b"TAG"), 1),
            (MismatchTriplet::new(*b"TAG", *b"TCG"), 1),
            (MismatchTriplet::new(*b"TAC", *b"TCC"), 1),
        ])
    );
    assert_eq!(
        count.insertion_length,
        HashMap::from([(1, 3), (2, 2), (3, 2), (4, 1)])
    );
    assert_eq!(
        count.deletion_length,
        HashMap::from([(1, 4), (2, 2), (3, 2), (4, 1)])
    );

    bam.read(&mut record).unwrap()?;
    assert_eq!(record.qname(), b"READ9");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;

    assert!(bam.read(&mut record).is_none());

    assert_eq!(count.total_sequenced_len, 1000);
    assert_eq!(*count.total_reference_triplet.get(b"NGA").unwrap(), 1);
    assert_eq!(*count.total_reference_triplet.get(b"TGN").unwrap(), 1);

    Ok(())
}

#[test]
fn test_add_bam() -> anyhow::Result<()> {
    let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
    let regions = Regions::create_from_fasta(&mut reference_fasta);
    let mut processor = SequencingErrorProcessor::new(10, None, reference_fasta, regions);

    let mut bam = bam::Reader::from_path("./testdata/demo1.bam")?;
    processor.add_bam(&mut bam, 0)?;
    let count = &processor.count;
    assert_eq!(
        count.mismatch,
        HashMap::from([
            (Mismatch::new(b'C', b'G'), 2),
            (Mismatch::new(b'T', b'A'), 2),
            (Mismatch::new(b'G', b'T'), 1),
            (Mismatch::new(b'T', b'G'), 3),
            (Mismatch::new(b'A', b'G'), 2),
            (Mismatch::new(b'A', b'C'), 3),
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
            (MismatchTriplet::new(*b"TTG", *b"TAG"), 1),
            (MismatchTriplet::new(*b"TAG", *b"TCG"), 1),
            (MismatchTriplet::new(*b"TAC", *b"TCC"), 1),
        ])
    );
    assert_eq!(
        count.insertion_length,
        HashMap::from([(1, 3), (2, 2), (3, 2), (4, 1)])
    );
    assert_eq!(
        count.deletion_length,
        HashMap::from([(1, 4), (2, 2), (3, 2), (4, 1)])
    );
    assert_eq!(count.softclip_length, HashMap::from([(9, 1), (6, 1)]));
    Ok(())
}

#[test]
fn test_add_cram() -> anyhow::Result<()> {
    let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
    let regions = Regions::create_from_fasta(&mut reference_fasta);
    let mut processor = SequencingErrorProcessor::new(10, None, reference_fasta, regions);

    let mut bam = bam::Reader::from_path("./testdata/demo1.cram")?;
    bam.set_reference("./testdata/ref/MT.fa")?;
    processor.add_bam(&mut bam, 0)?;
    let count = &processor.count;
    assert_eq!(
        count.mismatch,
        HashMap::from([
            (Mismatch::new(b'C', b'G'), 2),
            (Mismatch::new(b'T', b'A'), 2),
            (Mismatch::new(b'G', b'T'), 1),
            (Mismatch::new(b'T', b'G'), 3),
            (Mismatch::new(b'A', b'G'), 2),
            (Mismatch::new(b'A', b'C'), 3),
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
            (MismatchTriplet::new(*b"TTG", *b"TAG"), 1),
            (MismatchTriplet::new(*b"TAG", *b"TCG"), 1),
            (MismatchTriplet::new(*b"TAC", *b"TCC"), 1),
        ])
    );
    assert_eq!(
        count.insertion_length,
        HashMap::from([(1, 3), (2, 2), (3, 2), (4, 1)])
    );
    assert_eq!(
        count.deletion_length,
        HashMap::from([(1, 4), (2, 2), (3, 2), (4, 1)])
    );
    assert_eq!(count.softclip_length, HashMap::from([(9, 1), (6, 1)]));
    Ok(())
}

#[test]
fn test_add_sam() -> anyhow::Result<()> {
    let reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
    let regions = Regions::load_from_bed(File::open("testdata/target.bed")?)?;
    let mut processor = SequencingErrorProcessor::new(10, None, reference_fasta, regions);

    let mut bam = bam::Reader::from_path("./testdata/demo1.sam")?;
    bam.set_reference("./testdata/ref/MT.fa")?;
    processor.add_bam(&mut bam, 0)?;
    let count = &processor.count;
    assert_eq!(
        count.mismatch,
        HashMap::from([
            (Mismatch::new(b'C', b'G'), 1),
            (Mismatch::new(b'T', b'A'), 1),
        ])
    );
    assert_eq!(
        count.mismatch_triplet,
        HashMap::from([
            (MismatchTriplet::new(*b"CCC", *b"CGC"), 1),
            (MismatchTriplet::new(*b"TTT", *b"TAT"), 1),
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
    assert_eq!(count.softclip_length, HashMap::from([]));
    Ok(())
}

#[test]
fn test_add_record_with_known_variants() -> anyhow::Result<()> {
    let known_variants =
        rust_htslib::bcf::IndexedReader::from_path("./testdata/known_variants.vcf.gz")?;
    let mut reference_fasta = IndexedReader::from_file(&"./testdata/ref/MT.fa")?;
    let regions = Regions::create_from_fasta(&mut reference_fasta);
    let mut processor =
        SequencingErrorProcessor::new(10, Some(known_variants), reference_fasta, regions);

    let mut bam = bam::Reader::from_path("./testdata/demo1.bam")?;
    let mut record = bam::record::Record::new();

    bam.read(&mut record).unwrap()?;
    assert_eq!(record.qname(), b"READ8");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
    assert!(count.mismatch.is_empty());
    assert!(count.mismatch_triplet.is_empty());
    assert!(count.insertion_length.is_empty());
    assert!(count.deletion_length.is_empty());
    assert!(count.softclip_length.is_empty());
    assert_eq!(count.total_reference_len, 100);
    assert_eq!(count.total_sequenced_len, 100);

    bam.read(&mut record).unwrap()?;
    assert_eq!(record.qname(), b"READ1");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
    assert!(count.mismatch.is_empty());
    assert!(count.mismatch_triplet.is_empty());
    assert!(count.insertion_length.is_empty());
    assert!(count.deletion_length.is_empty());
    assert!(count.softclip_length.is_empty());
    assert_eq!(count.total_reference_len, 200);
    assert_eq!(count.total_sequenced_len, 200);

    bam.read(&mut record).unwrap()?;
    assert_eq!(record.qname(), b"READ6");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
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
    assert_eq!(record.qname(), b"READ10");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;
    assert_eq!(
        count.mismatch,
        HashMap::from([
            (Mismatch::new(b'C', b'G'), 2),
            (Mismatch::new(b'G', b'T'), 1),
            (Mismatch::new(b'T', b'G'), 3),
            (Mismatch::new(b'A', b'G'), 2),
            (Mismatch::new(b'A', b'C'), 3),
            (Mismatch::new(b'T', b'A'), 1),
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
            (MismatchTriplet::new(*b"TTG", *b"TAG"), 1),
            (MismatchTriplet::new(*b"TAG", *b"TCG"), 1),
            (MismatchTriplet::new(*b"TAC", *b"TCC"), 1),
        ])
    );
    assert_eq!(
        count.insertion_length,
        HashMap::from([(1, 3), (2, 1), (3, 2), (4, 1)])
    );
    assert_eq!(
        count.deletion_length,
        HashMap::from([(1, 3), (2, 2), (3, 2), (4, 1)])
    );

    bam.read(&mut record).unwrap()?;
    assert_eq!(record.qname(), b"READ9");
    processor.add_record(bam.header(), &record)?;
    let count = &processor.count;

    assert!(bam.read(&mut record).is_none());

    assert_eq!(count.total_sequenced_len, 1000);
    assert_eq!(*count.total_reference_triplet.get(b"NGA").unwrap(), 1);
    assert_eq!(*count.total_reference_triplet.get(b"TGN").unwrap(), 1);

    Ok(())
}

#[test]
fn test_regions() -> anyhow::Result<()> {
    let regions = Regions::load_from_bed(File::open("./testdata/target.bed")?)?;
    assert!(!regions.contains(b"MT", 9));
    assert!(regions.contains(b"MT", 10));
    assert!(regions.contains(b"MT", 1000));
    assert!(regions.contains(b"MT", 4999));
    assert!(!regions.contains(b"MT", 5000));
    assert!(!regions.contains(b"MT", 5001));

    Ok(())
}

#[test]
fn test_regions2() -> anyhow::Result<()> {
    let regions = Regions::load_from_bed(File::open("./testdata/region-test.bed")?)?;
    assert!(!regions.contains(b"chr1", 9));
    assert!(regions.contains(b"chr1", 10));
    assert!(regions.contains(b"chr1", 4999));
    assert!(!regions.contains(b"chr1", 5000));
    assert!(regions.contains(b"chr1", 7000));
    assert!(regions.contains(b"chr1", 9999));
    assert!(!regions.contains(b"chr1", 10000));
    assert!(regions.contains(b"chr1", 30000));
    assert!(regions.contains(b"chr1", 49999));
    assert!(!regions.contains(b"chr1", 50000));
    assert!(!regions.contains(b"chr2", 5999));
    assert!(regions.contains(b"chr2", 6000));
    assert!(regions.contains(b"chr2", 59999));
    assert!(!regions.contains(b"chr2", 60000));

    assert_eq!(
        regions,
        Regions {
            regions: HashMap::from([
                (
                    b"chr1".to_vec(),
                    vec![(10, 5000), (7000, 10000), (20000, 50000)]
                ),
                (b"chr2".to_vec(), vec![(6000, 60000)]),
                (b"chr3".to_vec(), vec![(10, 80)])
            ])
        }
    );

    Ok(())
}
