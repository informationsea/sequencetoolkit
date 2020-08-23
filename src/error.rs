use failure::{Backtrace, Context, Fail};
use std::fmt::{self, Display};
use std::num;

use super::geneannot::GeneAnnotError;
use super::vcfutils::error::VCFUtilsError;

#[derive(Debug)]
pub struct SequenceToolkitError {
    inner: failure::Context<SequenceToolkitErrorKind>,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Fail)]
pub enum SequenceToolkitErrorKind {
    #[fail(display = "VCFUtilsError")]
    VCFUtilsError,
    #[fail(display = "CSV Error")]
    CsvError,
    #[fail(display = "Gene Annotation Error")]
    GeneAnnotError,
    #[fail(display = "I/O Error")]
    IoError,
    #[fail(display = "VCFParse Error")]
    VCFError,
    #[fail(display = "Utf8 Error")]
    Utf8Error,
    #[fail(display = "Parse Int Error")]
    ParseIntError,
    #[fail(display = "Parse Float Error")]
    ParseFloatError,
    #[fail(display = "Error: {}", _0)]
    OtherError(&'static str),
    #[fail(display = "BED Parse Error: {}", _0)]
    BedParseError(u64),
}

impl Fail for SequenceToolkitError {
    fn cause(&self) -> Option<&dyn Fail> {
        self.inner.cause()
    }

    fn backtrace(&self) -> Option<&Backtrace> {
        self.inner.backtrace()
    }
}

impl Display for SequenceToolkitError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.inner, f)
    }
}

impl SequenceToolkitError {
    pub fn kind(&self) -> SequenceToolkitErrorKind {
        *self.inner.get_context()
    }
}

impl From<SequenceToolkitErrorKind> for SequenceToolkitError {
    fn from(kind: SequenceToolkitErrorKind) -> SequenceToolkitError {
        SequenceToolkitError {
            inner: Context::new(kind),
        }
    }
}

impl From<Context<SequenceToolkitErrorKind>> for SequenceToolkitError {
    fn from(inner: Context<SequenceToolkitErrorKind>) -> SequenceToolkitError {
        SequenceToolkitError { inner }
    }
}

impl From<VCFUtilsError> for SequenceToolkitError {
    fn from(e: VCFUtilsError) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::VCFUtilsError).into()
    }
}

impl From<GeneAnnotError> for SequenceToolkitError {
    fn from(e: GeneAnnotError) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::GeneAnnotError).into()
    }
}

impl From<std::io::Error> for SequenceToolkitError {
    fn from(e: std::io::Error) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::IoError).into()
    }
}

impl From<vcf::VCFError> for SequenceToolkitError {
    fn from(e: vcf::VCFError) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::VCFError).into()
    }
}

impl From<csv::Error> for SequenceToolkitError {
    fn from(e: csv::Error) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::CsvError).into()
    }
}

impl From<std::str::Utf8Error> for SequenceToolkitError {
    fn from(e: std::str::Utf8Error) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::Utf8Error).into()
    }
}

impl From<num::ParseIntError> for SequenceToolkitError {
    fn from(e: num::ParseIntError) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::ParseIntError).into()
    }
}

impl From<num::ParseFloatError> for SequenceToolkitError {
    fn from(e: num::ParseFloatError) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::ParseFloatError).into()
    }
}
