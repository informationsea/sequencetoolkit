use thiserror::Error;

#[derive(Debug, Error)]
pub enum GeneAnnotError {
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
    #[error("Parse Int Error: {0}")]
    ParseIntError(#[from] std::num::ParseIntError),
    #[error("Parse Float Error: {0}")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("Convert Error: {0}")]
    ConvertError(#[from] std::num::TryFromIntError),
    #[error("Serialize Error: {0}")]
    SerializeError(#[from] serde_json::Error),
    #[error("Bincode Serialize Error: {0}")]
    BincodeSerializeError(#[from] bincode::Error),
    #[error("refGene parse error at line {0}: {1}")]
    RefGeneParseError(u64, Box<GeneAnnotError>),
    #[error("Error: {0}")]
    OtherError(&'static str),
    #[error("HGVS position parse error")]
    HgvsPositionParseError,
    #[error("Error: {0}")]
    AnyhowError(#[from] anyhow::Error),
}
