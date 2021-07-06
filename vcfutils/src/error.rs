use failure::{Backtrace, Context, Fail};
use sequencetoolkit_common::{SequenceToolkitError, SequenceToolkitErrorKind};
use std::fmt::{self, Display};
use std::num;

#[derive(Debug)]
pub struct VCFUtilsError {
    inner: failure::Context<VCFUtilsErrorKind>,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Fail)]
pub enum VCFUtilsErrorKind {
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
    #[fail(display = "Genotype Parse Error")]
    GenotypeParseError,
    #[fail(display = "Parse Int Error")]
    ParseIntError,
    #[fail(display = "Parse Float Error")]
    ParseFloatError,
    #[fail(display = "Xlsx handling Error")]
    XlsxError,
    #[fail(display = "Error: {}", _0)]
    OtherError(&'static str),
}

impl Fail for VCFUtilsError {
    fn cause(&self) -> Option<&dyn Fail> {
        self.inner.cause()
    }

    fn backtrace(&self) -> Option<&Backtrace> {
        self.inner.backtrace()
    }
}

impl Display for VCFUtilsError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.inner, f)
    }
}

impl VCFUtilsError {
    pub fn kind(&self) -> VCFUtilsErrorKind {
        *self.inner.get_context()
    }
}

impl From<VCFUtilsErrorKind> for VCFUtilsError {
    fn from(kind: VCFUtilsErrorKind) -> VCFUtilsError {
        VCFUtilsError {
            inner: Context::new(kind),
        }
    }
}

impl From<Context<VCFUtilsErrorKind>> for VCFUtilsError {
    fn from(inner: Context<VCFUtilsErrorKind>) -> VCFUtilsError {
        VCFUtilsError { inner }
    }
}

impl From<vcf::VCFError> for VCFUtilsError {
    fn from(e: vcf::VCFError) -> VCFUtilsError {
        e.context(VCFUtilsErrorKind::VCFError).into()
    }
}

impl From<std::io::Error> for VCFUtilsError {
    fn from(e: std::io::Error) -> VCFUtilsError {
        e.context(VCFUtilsErrorKind::IoError).into()
    }
}

impl From<std::str::Utf8Error> for VCFUtilsError {
    fn from(e: std::str::Utf8Error) -> VCFUtilsError {
        e.context(VCFUtilsErrorKind::Utf8Error).into()
    }
}

impl From<csv::Error> for VCFUtilsError {
    fn from(e: csv::Error) -> VCFUtilsError {
        e.context(VCFUtilsErrorKind::CsvError).into()
    }
}

impl From<num::ParseIntError> for VCFUtilsError {
    fn from(e: num::ParseIntError) -> VCFUtilsError {
        e.context(VCFUtilsErrorKind::ParseIntError).into()
    }
}

impl From<num::ParseFloatError> for VCFUtilsError {
    fn from(e: num::ParseFloatError) -> VCFUtilsError {
        e.context(VCFUtilsErrorKind::ParseFloatError).into()
    }
}

impl From<fmt::Error> for VCFUtilsError {
    fn from(e: fmt::Error) -> VCFUtilsError {
        e.context(VCFUtilsErrorKind::FormatError).into()
    }
}

impl From<xlsxwriter::XlsxError> for VCFUtilsError {
    fn from(e: xlsxwriter::XlsxError) -> VCFUtilsError {
        e.context(VCFUtilsErrorKind::XlsxError).into()
    }
}

impl From<VCFUtilsError> for SequenceToolkitError {
    fn from(e: VCFUtilsError) -> SequenceToolkitError {
        e.context(SequenceToolkitErrorKind::VCFUtilsError).into()
    }
}
