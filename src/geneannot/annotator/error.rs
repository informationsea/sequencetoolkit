use failure::{Backtrace, Context, Fail};
use std::fmt::{self, Display};
use std::num;

#[derive(Debug)]
pub struct GeneAnnotError {
    inner: failure::Context<GeneAnnotErrorKind>,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Fail)]
pub enum GeneAnnotErrorKind {
    #[fail(display = "VCFParse Error")]
    VCFError,
    #[fail(display = "I/O Error")]
    IoError,
    #[fail(display = "Format Error")]
    FormatError,
    #[fail(display = "Utf8 Error")]
    Utf8Error,
    #[fail(display = "CSV Error")]
    CsvError,
    #[fail(display = "Parse Int Error")]
    ParseIntError,
    #[fail(display = "Parse Float Error")]
    ParseFloatError,
    #[fail(display = "Convert Error")]
    ConvertError,
    #[fail(display = "Serialize Error")]
    SerializeError,
    #[fail(display = "refGene parse error at line {}", _0)]
    RefGeneParseError(u64),
    #[fail(display = "Error: {}", _0)]
    OtherError(&'static str),
    #[fail(display = "HGVS position parse error")]
    HgvsPositionParseError,
}

impl Fail for GeneAnnotError {
    fn cause(&self) -> Option<&dyn Fail> {
        self.inner.cause()
    }

    fn backtrace(&self) -> Option<&Backtrace> {
        self.inner.backtrace()
    }
}

impl Display for GeneAnnotError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.inner, f)
    }
}

impl GeneAnnotError {
    pub fn kind(&self) -> GeneAnnotErrorKind {
        *self.inner.get_context()
    }
}

impl From<GeneAnnotErrorKind> for GeneAnnotError {
    fn from(kind: GeneAnnotErrorKind) -> GeneAnnotError {
        GeneAnnotError {
            inner: Context::new(kind),
        }
    }
}

impl From<Context<GeneAnnotErrorKind>> for GeneAnnotError {
    fn from(inner: Context<GeneAnnotErrorKind>) -> GeneAnnotError {
        GeneAnnotError { inner }
    }
}

impl From<vcf::VCFError> for GeneAnnotError {
    fn from(e: vcf::VCFError) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::VCFError).into()
    }
}

impl From<std::io::Error> for GeneAnnotError {
    fn from(e: std::io::Error) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::IoError).into()
    }
}

impl From<std::str::Utf8Error> for GeneAnnotError {
    fn from(e: std::str::Utf8Error) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::Utf8Error).into()
    }
}

impl From<csv::Error> for GeneAnnotError {
    fn from(e: csv::Error) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::CsvError).into()
    }
}

impl From<num::ParseIntError> for GeneAnnotError {
    fn from(e: num::ParseIntError) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::ParseIntError).into()
    }
}

impl From<num::ParseFloatError> for GeneAnnotError {
    fn from(e: num::ParseFloatError) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::ParseFloatError).into()
    }
}

impl From<fmt::Error> for GeneAnnotError {
    fn from(e: fmt::Error) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::FormatError).into()
    }
}

impl From<num::TryFromIntError> for GeneAnnotError {
    fn from(e: num::TryFromIntError) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::ConvertError).into()
    }
}

impl From<serde_json::Error> for GeneAnnotError {
    fn from(e: serde_json::Error) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::SerializeError).into()
    }
}

impl From<bincode::Error> for GeneAnnotError {
    fn from(e: bincode::Error) -> GeneAnnotError {
        e.context(GeneAnnotErrorKind::SerializeError).into()
    }
}
