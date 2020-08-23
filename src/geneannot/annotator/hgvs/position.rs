use crate::geneannot::annotator::models::{CdsPosition, TranscriptPosition};
use crate::geneannot::{GeneAnnotError, GeneAnnotErrorKind};

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::digit1;
use nom::combinator::{map, map_res, opt};
use nom::error::ErrorKind;
use nom::sequence::tuple;
use std::fmt;
use std::ops::{Deref, RangeInclusive};

#[derive(Debug, PartialEq, Eq, Copy, Clone, PartialOrd, Ord)]
pub enum ParsedPosition {
    CdsPosition(CdsPosition),
    TranscriptPosition(TranscriptPosition),
    GenomePosition(u64),
}

fn write_partial_transcript_position(
    formatter: &mut fmt::Formatter<'_>,
    transcript_position: TranscriptPosition,
) -> fmt::Result {
    match transcript_position {
        TranscriptPosition::Exon(x) => write!(formatter, "{}", x + 1),
        TranscriptPosition::BeforeExon(x, o) => write!(formatter, "{}-{}", x + 1, o),
        TranscriptPosition::AfterExon(x, o) => write!(formatter, "{}+{}", x + 1, o),
    }
}

impl fmt::Display for CdsPosition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "c.")?;
        write_partial_cds_position(f, *self)
    }
}

pub struct TranscriptPositionRangeInclusive(RangeInclusive<TranscriptPosition>);

impl Deref for TranscriptPositionRangeInclusive {
    type Target = RangeInclusive<TranscriptPosition>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for TranscriptPositionRangeInclusive {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "n.")?;
        write_partial_transcript_position(f, *self.start())?;
        write!(f, "_")?;
        write_partial_transcript_position(f, *self.end())
    }
}

impl Into<TranscriptPositionRangeInclusive> for RangeInclusive<TranscriptPosition> {
    fn into(self) -> TranscriptPositionRangeInclusive {
        TranscriptPositionRangeInclusive(self)
    }
}

impl fmt::Display for TranscriptPosition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "n.")?;
        write_partial_transcript_position(f, *self)
    }
}

fn write_partial_cds_position(
    formatter: &mut fmt::Formatter<'_>,
    cds_position: CdsPosition,
) -> fmt::Result {
    match cds_position {
        CdsPosition::Cds(x) => write_partial_transcript_position(formatter, x),
        CdsPosition::BeforeCds(x) => {
            write!(formatter, "-")?;
            write_partial_transcript_position(formatter, x)
        }
        CdsPosition::AfterCds(x) => {
            write!(formatter, "*")?;
            write_partial_transcript_position(formatter, x)
        }
    }
}

fn partial_transcript_position_parser(
    text: &str,
) -> Result<(&str, TranscriptPosition), nom::Err<(&str, ErrorKind)>> {
    map_res::<_, _, _, (&str, ErrorKind), GeneAnnotError, _, _>(
        tuple((digit1, opt(tuple((alt((tag("+"), tag("-"))), digit1))))),
        |r| {
            let base_position = r.0.parse::<u64>()? - 1;
            match r.1 {
                Some(("+", x)) => Ok(TranscriptPosition::AfterExon(base_position, x.parse()?)),
                Some(("-", x)) => Ok(TranscriptPosition::BeforeExon(base_position, x.parse()?)),
                None => Ok(TranscriptPosition::Exon(base_position)),
                _ => Err(GeneAnnotErrorKind::HgvsPositionParseError.into()),
            }
        },
    )(text)
}

pub fn parse_hgvs_position(text: &str) -> Result<ParsedPosition, GeneAnnotError> {
    match alt((
        map_res::<_, _, _, (&str, ErrorKind), GeneAnnotError, _, _>(
            tuple((tag("g."), digit1)),
            |r| Ok(ParsedPosition::GenomePosition(r.1.parse::<u64>()? - 1)),
        ),
        map::<_, _, _, (&str, ErrorKind), _, _>(
            tuple((tag("n."), partial_transcript_position_parser)),
            |r| ParsedPosition::TranscriptPosition(r.1),
        ),
        map_res::<_, _, _, (&str, ErrorKind), GeneAnnotError, _, _>(
            tuple((
                tag("c."),
                opt(alt((tag("*"), tag("-")))),
                partial_transcript_position_parser,
            )),
            |r| match r.1 {
                Some("*") => Ok(ParsedPosition::CdsPosition(CdsPosition::AfterCds(r.2))),
                Some("-") => Ok(ParsedPosition::CdsPosition(CdsPosition::BeforeCds(r.2))),
                None => Ok(ParsedPosition::CdsPosition(CdsPosition::Cds(r.2))),
                _ => Err(GeneAnnotErrorKind::HgvsPositionParseError.into()),
            },
        ),
    ))(text)
    {
        Ok(("", result)) => Ok(result),
        _ => Err(GeneAnnotErrorKind::HgvsPositionParseError.into()),
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_write_hgvs_cds_position() {
        assert_eq!(
            format!("{}", CdsPosition::Cds(TranscriptPosition::Exon(10))),
            "c.11"
        );
        assert_eq!(
            format!("{}", CdsPosition::BeforeCds(TranscriptPosition::Exon(10))),
            "c.-11"
        );
        assert_eq!(
            format!("{}", CdsPosition::AfterCds(TranscriptPosition::Exon(10))),
            "c.*11"
        );
        assert_eq!(
            format!(
                "{}",
                CdsPosition::Cds(TranscriptPosition::BeforeExon(10, 2))
            ),
            "c.11-2"
        );
        assert_eq!(
            format!(
                "{}",
                CdsPosition::BeforeCds(TranscriptPosition::BeforeExon(10, 2))
            ),
            "c.-11-2"
        );
        assert_eq!(
            format!(
                "{}",
                CdsPosition::AfterCds(TranscriptPosition::BeforeExon(10, 2))
            ),
            "c.*11-2"
        );
        assert_eq!(
            format!("{}", CdsPosition::Cds(TranscriptPosition::AfterExon(10, 3))),
            "c.11+3"
        );
        assert_eq!(
            format!(
                "{}",
                CdsPosition::BeforeCds(TranscriptPosition::AfterExon(10, 3))
            ),
            "c.-11+3"
        );
        assert_eq!(
            format!(
                "{}",
                CdsPosition::AfterCds(TranscriptPosition::AfterExon(10, 3))
            ),
            "c.*11+3"
        );
    }

    #[test]
    fn test_parse_hgvs_position() {
        assert_eq!(
            parse_hgvs_position("g.123").unwrap(),
            ParsedPosition::GenomePosition(122)
        );

        assert_eq!(
            parse_hgvs_position("n.123").unwrap(),
            ParsedPosition::TranscriptPosition(TranscriptPosition::Exon(122))
        );
        assert_eq!(
            parse_hgvs_position("n.123+4").unwrap(),
            ParsedPosition::TranscriptPosition(TranscriptPosition::AfterExon(122, 4))
        );
        assert_eq!(
            parse_hgvs_position("n.123-5").unwrap(),
            ParsedPosition::TranscriptPosition(TranscriptPosition::BeforeExon(122, 5))
        );

        assert_eq!(
            parse_hgvs_position("c.123").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::Cds(TranscriptPosition::Exon(122)))
        );
        assert_eq!(
            parse_hgvs_position("c.123+4").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::Cds(TranscriptPosition::AfterExon(122, 4)))
        );
        assert_eq!(
            parse_hgvs_position("c.123-5").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::Cds(TranscriptPosition::BeforeExon(122, 5)))
        );

        assert_eq!(
            parse_hgvs_position("c.*123").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::AfterCds(TranscriptPosition::Exon(122)))
        );
        assert_eq!(
            parse_hgvs_position("c.*123+4").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::AfterCds(TranscriptPosition::AfterExon(
                122, 4
            )))
        );
        assert_eq!(
            parse_hgvs_position("c.*123-5").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::AfterCds(TranscriptPosition::BeforeExon(
                122, 5
            )))
        );

        assert_eq!(
            parse_hgvs_position("c.-123").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::BeforeCds(TranscriptPosition::Exon(122)))
        );
        assert_eq!(
            parse_hgvs_position("c.-123+4").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::BeforeCds(TranscriptPosition::AfterExon(
                122, 4
            )))
        );
        assert_eq!(
            parse_hgvs_position("c.-123-5").unwrap(),
            ParsedPosition::CdsPosition(CdsPosition::BeforeCds(TranscriptPosition::BeforeExon(
                122, 5
            )))
        );
    }
}
