use thiserror::Error;

#[derive(Debug, Error)]
pub enum VCFUtilsError {
    #[error("VCFParse Error: {0}")]
    VCFError(#[from] vcf::VCFError),
    #[error("I/O Error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Format Error: {0}")]
    FormatError(#[from] std::fmt::Error),
    #[error("Utf8 Error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("CSV Error: {0}")]
    CsvError(#[from] csv::Error),
    #[error("Genotype Parse Error")]
    GenotypeParseError,
    #[error("Parse Int Error: {0}")]
    ParseIntError(#[from] std::num::ParseIntError),
    #[error("Parse Float Error: {0}")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("Xlsx handling Error: {0}")]
    XlsxError(#[from] xlsxwriter::XlsxError),
    #[error("Error: {0}")]
    OtherError(&'static str),
}
