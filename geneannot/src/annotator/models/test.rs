use super::*;

#[test]
fn test_cds_position_ordering() {
    assert_eq!(
        CdsPosition::Cds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::Cds(TranscriptPosition::Exon(10))),
        Ordering::Equal
    );
    assert_eq!(
        CdsPosition::Cds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::Cds(TranscriptPosition::Exon(9))),
        Ordering::Greater
    );
    assert_eq!(
        CdsPosition::Cds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::Cds(TranscriptPosition::Exon(11))),
        Ordering::Less
    );
    assert_eq!(
        CdsPosition::Cds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::BeforeCds(TranscriptPosition::Exon(10))),
        Ordering::Greater
    );
    assert_eq!(
        CdsPosition::Cds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::AfterCds(TranscriptPosition::Exon(10))),
        Ordering::Less
    );

    assert_eq!(
        CdsPosition::BeforeCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::BeforeCds(TranscriptPosition::Exon(10))),
        Ordering::Equal
    );
    assert_eq!(
        CdsPosition::BeforeCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::BeforeCds(TranscriptPosition::Exon(9))),
        Ordering::Less
    );
    assert_eq!(
        CdsPosition::BeforeCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::BeforeCds(TranscriptPosition::Exon(11))),
        Ordering::Greater
    );
    assert_eq!(
        CdsPosition::BeforeCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::Cds(TranscriptPosition::Exon(10))),
        Ordering::Less
    );
    assert_eq!(
        CdsPosition::BeforeCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::AfterCds(TranscriptPosition::Exon(10))),
        Ordering::Less
    );

    assert_eq!(
        CdsPosition::AfterCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::AfterCds(TranscriptPosition::Exon(10))),
        Ordering::Equal
    );
    assert_eq!(
        CdsPosition::AfterCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::AfterCds(TranscriptPosition::Exon(11))),
        Ordering::Less
    );
    assert_eq!(
        CdsPosition::AfterCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::AfterCds(TranscriptPosition::Exon(9))),
        Ordering::Greater
    );
    assert_eq!(
        CdsPosition::AfterCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::Cds(TranscriptPosition::Exon(10))),
        Ordering::Greater
    );
    assert_eq!(
        CdsPosition::AfterCds(TranscriptPosition::Exon(10))
            .cmp(&CdsPosition::BeforeCds(TranscriptPosition::Exon(10))),
        Ordering::Greater
    );
}

#[test]
fn test_transcript_position_ordering() {
    assert_eq!(
        TranscriptPosition::Exon(10).cmp(&TranscriptPosition::BeforeExon(11, 1)),
        Ordering::Less
    );
    assert_eq!(
        TranscriptPosition::Exon(10).cmp(&TranscriptPosition::AfterExon(9, 1)),
        Ordering::Greater
    );
    assert_eq!(
        TranscriptPosition::AfterExon(10, 1).cmp(&TranscriptPosition::BeforeExon(11, 1)),
        Ordering::Less
    );
    assert_eq!(
        TranscriptPosition::AfterExon(10, 1).cmp(&TranscriptPosition::Exon(11)),
        Ordering::Less
    );
    assert_eq!(
        TranscriptPosition::BeforeExon(10, 1).cmp(&TranscriptPosition::AfterExon(9, 1)),
        Ordering::Greater
    );
    assert_eq!(
        TranscriptPosition::BeforeExon(10, 1).cmp(&TranscriptPosition::Exon(9)),
        Ordering::Greater
    );

    assert_eq!(
        TranscriptPosition::Exon(10).cmp(&TranscriptPosition::Exon(10)),
        Ordering::Equal
    );
    assert_eq!(
        TranscriptPosition::Exon(10).cmp(&TranscriptPosition::BeforeExon(10, 1)),
        Ordering::Greater
    );
    assert_eq!(
        TranscriptPosition::Exon(10).cmp(&TranscriptPosition::AfterExon(10, 1)),
        Ordering::Less
    );
    assert!(TranscriptPosition::Exon(10) < TranscriptPosition::AfterExon(10, 1));

    assert_eq!(
        TranscriptPosition::BeforeExon(10, 2).cmp(&TranscriptPosition::Exon(10)),
        Ordering::Less
    );
    assert_eq!(
        TranscriptPosition::BeforeExon(10, 2).cmp(&TranscriptPosition::BeforeExon(10, 1)),
        Ordering::Less
    );
    assert_eq!(
        TranscriptPosition::BeforeExon(10, 2).cmp(&TranscriptPosition::BeforeExon(10, 2)),
        Ordering::Equal
    );
    assert_eq!(
        TranscriptPosition::BeforeExon(10, 2).cmp(&TranscriptPosition::BeforeExon(10, 3)),
        Ordering::Greater
    );
    assert_eq!(
        TranscriptPosition::BeforeExon(10, 2).cmp(&TranscriptPosition::AfterExon(10, 1)),
        Ordering::Less
    );

    assert_eq!(
        TranscriptPosition::AfterExon(10, 2).cmp(&TranscriptPosition::Exon(10)),
        Ordering::Greater
    );
    assert_eq!(
        TranscriptPosition::AfterExon(10, 2).cmp(&TranscriptPosition::BeforeExon(10, 1)),
        Ordering::Greater
    );
    assert_eq!(
        TranscriptPosition::AfterExon(10, 2).cmp(&TranscriptPosition::AfterExon(10, 1)),
        Ordering::Greater
    );
    assert_eq!(
        TranscriptPosition::AfterExon(10, 2).cmp(&TranscriptPosition::AfterExon(10, 2)),
        Ordering::Equal
    );
    assert_eq!(
        TranscriptPosition::AfterExon(10, 2).cmp(&TranscriptPosition::AfterExon(10, 3)),
        Ordering::Less
    );
}

#[test]
fn test_exon_position_forward() {
    let exon = Exon {
        start: 10,
        end: 20,
        annotations: HashMap::new(),
    };

    assert_eq!(
        ExonPosition::BeforeExon(5),
        exon.exon_position(5, Strand::Forward)
    );
    assert_eq!(
        ExonPosition::BeforeExon(1),
        exon.exon_position(9, Strand::Forward)
    );
    assert_eq!(
        ExonPosition::InExon(0),
        exon.exon_position(10, Strand::Forward)
    );
    assert_eq!(
        ExonPosition::InExon(9),
        exon.exon_position(19, Strand::Forward)
    );
    assert_eq!(
        ExonPosition::AfterExon(1),
        exon.exon_position(20, Strand::Forward)
    );
    assert_eq!(
        ExonPosition::AfterExon(5),
        exon.exon_position(24, Strand::Forward)
    );
}

#[test]
fn test_exon_position_reverse() {
    let exon = Exon {
        start: 10,
        end: 20,
        annotations: HashMap::new(),
    };

    assert_eq!(
        ExonPosition::AfterExon(5),
        exon.exon_position(5, Strand::Reverse)
    );
    assert_eq!(
        ExonPosition::AfterExon(1),
        exon.exon_position(9, Strand::Reverse)
    );
    assert_eq!(
        ExonPosition::InExon(9),
        exon.exon_position(10, Strand::Reverse)
    );
    assert_eq!(
        ExonPosition::InExon(0),
        exon.exon_position(19, Strand::Reverse)
    );
    assert_eq!(
        ExonPosition::BeforeExon(1),
        exon.exon_position(20, Strand::Reverse)
    );
    assert_eq!(
        ExonPosition::BeforeExon(5),
        exon.exon_position(24, Strand::Reverse)
    );
}

#[test]
fn test_transcript_position_forward() {
    let transcript_forward = NoncodingTranscript {
        id: "foo".to_string(),
        chromosome_index: 0,
        strand: Strand::Forward,
        start: 100,
        end: 300,
        annotations: HashMap::new(),
        exons: vec![
            Exon {
                start: 100,
                end: 150,
                annotations: HashMap::new(),
            },
            Exon {
                start: 200,
                end: 230,
                annotations: HashMap::new(),
            },
            Exon {
                start: 250,
                end: 300,
                annotations: HashMap::new(),
            },
        ],
    };
    assert_eq!(
        transcript_forward.transcript_position(99),
        TranscriptPosition::BeforeExon(0, 1)
    );
    assert_eq!(
        transcript_forward.transcript_position(100),
        TranscriptPosition::Exon(0)
    );
    assert_eq!(
        transcript_forward.transcript_position(149),
        TranscriptPosition::Exon(49)
    );
    assert_eq!(
        transcript_forward.transcript_position(150),
        TranscriptPosition::AfterExon(49, 1)
    );
    assert_eq!(
        transcript_forward.transcript_position(199),
        TranscriptPosition::BeforeExon(50, 1)
    );
    assert_eq!(
        transcript_forward.transcript_position(200),
        TranscriptPosition::Exon(50)
    );
    assert_eq!(
        transcript_forward.transcript_position(300),
        TranscriptPosition::AfterExon(129, 1)
    );
}

#[test]
fn test_transcript_position_reverse() {
    let transcript_forward = NoncodingTranscript {
        id: "foo".to_string(),
        chromosome_index: 0,
        strand: Strand::Reverse,
        start: 100,
        end: 300,
        annotations: HashMap::new(),
        exons: vec![
            Exon {
                start: 100,
                end: 150,
                annotations: HashMap::new(),
            },
            Exon {
                start: 200,
                end: 230,
                annotations: HashMap::new(),
            },
            Exon {
                start: 250,
                end: 300,
                annotations: HashMap::new(),
            },
        ],
    };
    assert_eq!(
        transcript_forward.transcript_position(99),
        TranscriptPosition::AfterExon(129, 1)
    );
    assert_eq!(
        transcript_forward.transcript_position(100),
        TranscriptPosition::Exon(129)
    );
    assert_eq!(
        transcript_forward.transcript_position(149),
        TranscriptPosition::Exon(80)
    );
    assert_eq!(
        transcript_forward.transcript_position(150),
        TranscriptPosition::BeforeExon(80, 1)
    );
    assert_eq!(
        transcript_forward.transcript_position(199),
        TranscriptPosition::AfterExon(79, 1)
    );
    assert_eq!(
        transcript_forward.transcript_position(200),
        TranscriptPosition::Exon(79)
    );
    assert_eq!(
        transcript_forward.transcript_position(229),
        TranscriptPosition::Exon(50)
    );
    assert_eq!(
        transcript_forward.transcript_position(300),
        TranscriptPosition::BeforeExon(0, 1)
    );
}

#[test]
fn test_transcript_genome_position_forward() {
    let transcript_forward = NoncodingTranscript {
        id: "foo".to_string(),
        chromosome_index: 0,
        strand: Strand::Forward,
        start: 100,
        end: 300,
        annotations: HashMap::new(),
        exons: vec![
            Exon {
                start: 100,
                end: 150,
                annotations: HashMap::new(),
            },
            Exon {
                start: 200,
                end: 230,
                annotations: HashMap::new(),
            },
            Exon {
                start: 250,
                end: 300,
                annotations: HashMap::new(),
            },
        ],
    };
    assert_eq!(
        99,
        transcript_forward.genome_position(TranscriptPosition::BeforeExon(0, 1))
    );
    assert_eq!(
        100,
        transcript_forward.genome_position(TranscriptPosition::Exon(0))
    );
    assert_eq!(
        149,
        transcript_forward.genome_position(TranscriptPosition::Exon(49))
    );
    assert_eq!(
        150,
        transcript_forward.genome_position(TranscriptPosition::AfterExon(49, 1))
    );
    assert_eq!(
        199,
        transcript_forward.genome_position(TranscriptPosition::BeforeExon(50, 1))
    );
    assert_eq!(
        200,
        transcript_forward.genome_position(TranscriptPosition::Exon(50))
    );
    assert_eq!(
        300,
        transcript_forward.genome_position(TranscriptPosition::AfterExon(129, 1))
    );
}

#[test]
fn test_genome_position_reverse() {
    let transcript_forward = NoncodingTranscript {
        id: "foo".to_string(),
        chromosome_index: 0,
        strand: Strand::Reverse,
        start: 100,
        end: 300,
        annotations: HashMap::new(),
        exons: vec![
            Exon {
                start: 100,
                end: 150,
                annotations: HashMap::new(),
            },
            Exon {
                start: 200,
                end: 230,
                annotations: HashMap::new(),
            },
            Exon {
                start: 250,
                end: 300,
                annotations: HashMap::new(),
            },
        ],
    };
    assert_eq!(
        99,
        transcript_forward.genome_position(TranscriptPosition::AfterExon(129, 1))
    );
    assert_eq!(
        100,
        transcript_forward.genome_position(TranscriptPosition::Exon(129))
    );
    assert_eq!(
        149,
        transcript_forward.genome_position(TranscriptPosition::Exon(80))
    );
    assert_eq!(
        150,
        transcript_forward.genome_position(TranscriptPosition::BeforeExon(80, 1))
    );
    assert_eq!(
        199,
        transcript_forward.genome_position(TranscriptPosition::AfterExon(79, 1))
    );
    assert_eq!(
        200,
        transcript_forward.genome_position(TranscriptPosition::Exon(79))
    );
    assert_eq!(
        229,
        transcript_forward.genome_position(TranscriptPosition::Exon(50))
    );
    assert_eq!(
        300,
        transcript_forward.genome_position(TranscriptPosition::BeforeExon(0, 1))
    );
}

fn create_transcript_forward() -> CodingTranscript {
    CodingTranscript {
        id: "foo".to_string(),
        chromosome_index: 0,
        strand: Strand::Forward,
        start: 100,
        end: 300,
        cds_start: 210,
        cds_end: 270,
        cds_start_status: CdsStatus::Complete,
        cds_end_status: CdsStatus::Complete,
        annotations: HashMap::new(),
        exons: vec![
            Exon {
                start: 100,
                end: 150,
                annotations: HashMap::new(),
            },
            Exon {
                start: 200,
                end: 230,
                annotations: HashMap::new(),
            },
            Exon {
                start: 250,
                end: 300,
                annotations: HashMap::new(),
            },
            Exon {
                start: 350,
                end: 400,
                annotations: HashMap::new(),
            },
        ],
    }
}

fn create_transcript_reverse() -> CodingTranscript {
    CodingTranscript {
        id: "foo".to_string(),
        chromosome_index: 0,
        strand: Strand::Reverse,
        start: 100,
        end: 300,
        cds_start: 210,
        cds_end: 270,
        cds_start_status: CdsStatus::Complete,
        cds_end_status: CdsStatus::Complete,
        annotations: HashMap::new(),
        exons: vec![
            Exon {
                start: 100,
                end: 150,
                annotations: HashMap::new(),
            },
            Exon {
                start: 200,
                end: 230,
                annotations: HashMap::new(),
            },
            Exon {
                start: 250,
                end: 300,
                annotations: HashMap::new(),
            },
            Exon {
                start: 350,
                end: 400,
                annotations: HashMap::new(),
            },
        ],
    }
}

#[test]
fn test_cds_offset() {
    let transcript_forward = create_transcript_forward();
    assert_eq!(transcript_forward.cds_offset(), Some(60));

    let transcript_reverse = create_transcript_reverse();
    assert_eq!(transcript_reverse.cds_offset(), Some(80));
}

#[test]
fn test_len() {
    let transcript_forward = create_transcript_forward();
    assert_eq!(transcript_forward.len(), 50 + 30 + 50 + 50);
    assert_eq!(transcript_forward.cds_len(), Some(20 + 20));

    let transcript_reverse = create_transcript_reverse();
    assert_eq!(transcript_reverse.len(), 50 + 30 + 50 + 50);
    assert_eq!(transcript_reverse.cds_len(), Some(20 + 20));
}

#[test]
fn test_cds_position_forward() {
    let transcript_forward = create_transcript_forward();
    assert_eq!(
        transcript_forward.cds_position(99),
        Some(CdsPosition::BeforeCds(TranscriptPosition::BeforeExon(
            59, 1
        )))
    );
    assert_eq!(
        transcript_forward.cds_position(100),
        Some(CdsPosition::BeforeCds(TranscriptPosition::Exon(59)))
    );
    assert_eq!(
        transcript_forward.cds_position(149),
        Some(CdsPosition::BeforeCds(TranscriptPosition::Exon(10)))
    );
    assert_eq!(
        transcript_forward.cds_position(150),
        Some(CdsPosition::BeforeCds(TranscriptPosition::AfterExon(10, 1)))
    );
    assert_eq!(
        transcript_forward.cds_position(199),
        Some(CdsPosition::BeforeCds(TranscriptPosition::BeforeExon(9, 1)))
    );
    assert_eq!(
        transcript_forward.cds_position(200),
        Some(CdsPosition::BeforeCds(TranscriptPosition::Exon(9)))
    );
    assert_eq!(
        transcript_forward.cds_position(209),
        Some(CdsPosition::BeforeCds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        transcript_forward.cds_position(210),
        Some(CdsPosition::Cds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        transcript_forward.cds_position(229),
        Some(CdsPosition::Cds(TranscriptPosition::Exon(19)))
    );
    assert_eq!(
        transcript_forward.cds_position(230),
        Some(CdsPosition::Cds(TranscriptPosition::AfterExon(19, 1)))
    );
    assert_eq!(
        transcript_forward.cds_position(249),
        Some(CdsPosition::Cds(TranscriptPosition::BeforeExon(20, 1)))
    );
    assert_eq!(
        transcript_forward.cds_position(250),
        Some(CdsPosition::Cds(TranscriptPosition::Exon(20)))
    );
    assert_eq!(
        transcript_forward.cds_position(269),
        Some(CdsPosition::Cds(TranscriptPosition::Exon(39)))
    );
    assert_eq!(
        transcript_forward.cds_position(270),
        Some(CdsPosition::AfterCds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        transcript_forward.cds_position(299),
        Some(CdsPosition::AfterCds(TranscriptPosition::Exon(29)))
    );
    assert_eq!(
        transcript_forward.cds_position(300),
        Some(CdsPosition::AfterCds(TranscriptPosition::AfterExon(29, 1)))
    );
    assert_eq!(
        transcript_forward.cds_position(349),
        Some(CdsPosition::AfterCds(TranscriptPosition::BeforeExon(30, 1)))
    );
    assert_eq!(
        transcript_forward.cds_position(350),
        Some(CdsPosition::AfterCds(TranscriptPosition::Exon(30)))
    );
    assert_eq!(
        transcript_forward.cds_position(399),
        Some(CdsPosition::AfterCds(TranscriptPosition::Exon(79)))
    );
    assert_eq!(
        transcript_forward.cds_position(400),
        Some(CdsPosition::AfterCds(TranscriptPosition::AfterExon(79, 1)))
    );
}

#[test]
fn test_cds_position_reverse() {
    let transcript_reverse = create_transcript_reverse();
    assert_eq!(
        transcript_reverse.cds_position(99),
        Some(CdsPosition::AfterCds(TranscriptPosition::AfterExon(59, 1)))
    );
    assert_eq!(
        transcript_reverse.cds_position(100),
        Some(CdsPosition::AfterCds(TranscriptPosition::Exon(59)))
    );
    assert_eq!(
        transcript_reverse.cds_position(149),
        Some(CdsPosition::AfterCds(TranscriptPosition::Exon(10)))
    );
    assert_eq!(
        transcript_reverse.cds_position(150),
        Some(CdsPosition::AfterCds(TranscriptPosition::BeforeExon(10, 1)))
    );
    assert_eq!(
        transcript_reverse.cds_position(199),
        Some(CdsPosition::AfterCds(TranscriptPosition::AfterExon(9, 1)))
    );
    assert_eq!(
        transcript_reverse.cds_position(200),
        Some(CdsPosition::AfterCds(TranscriptPosition::Exon(9)))
    );
    assert_eq!(
        transcript_reverse.cds_position(209),
        Some(CdsPosition::AfterCds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        transcript_reverse.cds_position(210),
        Some(CdsPosition::Cds(TranscriptPosition::Exon(39)))
    );
    assert_eq!(
        transcript_reverse.cds_position(229),
        Some(CdsPosition::Cds(TranscriptPosition::Exon(20)))
    );
    assert_eq!(
        transcript_reverse.cds_position(230),
        Some(CdsPosition::Cds(TranscriptPosition::BeforeExon(20, 1)))
    );
    assert_eq!(
        transcript_reverse.cds_position(249),
        Some(CdsPosition::Cds(TranscriptPosition::AfterExon(19, 1)))
    );
    assert_eq!(
        transcript_reverse.cds_position(250),
        Some(CdsPosition::Cds(TranscriptPosition::Exon(19)))
    );
    assert_eq!(
        transcript_reverse.cds_position(269),
        Some(CdsPosition::Cds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        transcript_reverse.cds_position(270),
        Some(CdsPosition::BeforeCds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        transcript_reverse.cds_position(299),
        Some(CdsPosition::BeforeCds(TranscriptPosition::Exon(29)))
    );
    assert_eq!(
        transcript_reverse.cds_position(300),
        Some(CdsPosition::BeforeCds(TranscriptPosition::BeforeExon(
            29, 1
        )))
    );
    assert_eq!(
        transcript_reverse.cds_position(349),
        Some(CdsPosition::BeforeCds(TranscriptPosition::AfterExon(30, 1)))
    );
    assert_eq!(
        transcript_reverse.cds_position(350),
        Some(CdsPosition::BeforeCds(TranscriptPosition::Exon(30)))
    );
    assert_eq!(
        transcript_reverse.cds_position(399),
        Some(CdsPosition::BeforeCds(TranscriptPosition::Exon(79)))
    );
    assert_eq!(
        transcript_reverse.cds_position(400),
        Some(CdsPosition::BeforeCds(TranscriptPosition::BeforeExon(
            79, 1
        )))
    );
}

#[test]
fn test_genome_position_from_cds_forward() {
    let transcript_forward = create_transcript_forward();
    assert_eq!(
        transcript_forward.cds_position(99),
        Some(CdsPosition::BeforeCds(TranscriptPosition::BeforeExon(
            59, 1
        )))
    );
    assert_eq!(
        Some(100),
        transcript_forward
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::Exon(59)))
    );
    assert_eq!(
        Some(149),
        transcript_forward
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::Exon(10)))
    );
    assert_eq!(
        Some(150),
        transcript_forward
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::AfterExon(10, 1)))
    );
    assert_eq!(
        Some(199),
        transcript_forward
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::BeforeExon(9, 1)))
    );
    assert_eq!(
        Some(200),
        transcript_forward
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::Exon(9)))
    );
    assert_eq!(
        Some(209),
        transcript_forward
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        Some(210),
        transcript_forward.genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        Some(229),
        transcript_forward.genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::Exon(19)))
    );
    assert_eq!(
        Some(230),
        transcript_forward
            .genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::AfterExon(19, 1)))
    );
    assert_eq!(
        Some(249),
        transcript_forward
            .genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::BeforeExon(20, 1)))
    );
    assert_eq!(
        Some(250),
        transcript_forward.genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::Exon(20)))
    );
    assert_eq!(
        Some(269),
        transcript_forward.genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::Exon(39)))
    );
    assert_eq!(
        Some(270),
        transcript_forward
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        Some(299),
        transcript_forward
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::Exon(29)))
    );
    assert_eq!(
        Some(300),
        transcript_forward
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::AfterExon(29, 1)))
    );
    assert_eq!(
        Some(349),
        transcript_forward
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::BeforeExon(30, 1)))
    );
    assert_eq!(
        Some(350),
        transcript_forward
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::Exon(30)))
    );
    assert_eq!(
        Some(399),
        transcript_forward
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::Exon(79)))
    );
    assert_eq!(
        Some(400),
        transcript_forward
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::AfterExon(79, 1)))
    );
}

#[test]
fn test_genome_position_from_cds_reverse() {
    let transcript_reverse = create_transcript_reverse();
    assert_eq!(
        Some(99),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::AfterExon(59, 1)))
    );
    assert_eq!(
        Some(100),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::Exon(59)))
    );
    assert_eq!(
        Some(149),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::Exon(10)))
    );
    assert_eq!(
        Some(150),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::BeforeExon(10, 1)))
    );
    assert_eq!(
        Some(199),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::AfterExon(9, 1)))
    );
    assert_eq!(
        Some(200),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::Exon(9)))
    );
    assert_eq!(
        Some(209),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::AfterCds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        Some(210),
        transcript_reverse.genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::Exon(39)))
    );
    assert_eq!(
        Some(229),
        transcript_reverse.genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::Exon(20)))
    );
    assert_eq!(
        Some(230),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::BeforeExon(20, 1)))
    );
    assert_eq!(
        Some(249),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::AfterExon(19, 1)))
    );
    assert_eq!(
        Some(250),
        transcript_reverse.genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::Exon(19)))
    );
    assert_eq!(
        Some(269),
        transcript_reverse.genome_position_from_cds(CdsPosition::Cds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        Some(270),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::Exon(0)))
    );
    assert_eq!(
        Some(299),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::Exon(29)))
    );
    assert_eq!(
        Some(300),
        transcript_reverse.genome_position_from_cds(CdsPosition::BeforeCds(
            TranscriptPosition::BeforeExon(29, 1)
        ))
    );
    assert_eq!(
        Some(349),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::AfterExon(30, 1)))
    );
    assert_eq!(
        Some(350),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::Exon(30)))
    );
    assert_eq!(
        Some(399),
        transcript_reverse
            .genome_position_from_cds(CdsPosition::BeforeCds(TranscriptPosition::Exon(79)))
    );
    assert_eq!(
        Some(400),
        transcript_reverse.genome_position_from_cds(CdsPosition::BeforeCds(
            TranscriptPosition::BeforeExon(79, 1)
        ))
    );
}
