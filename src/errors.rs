use thiserror::Error;

// define Error types
#[derive(Debug, Error)]
pub enum FqkitError {
    #[error("Stdin not detected")]
    StdinNotDetected,

    #[error("Failed to open file: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Invalid output dir: {0}")]
    InvalidOutputDir(String),

    #[error("EmptyQualRecord")]
    EmptyQualRecord,

    #[error("ThreadPoolBuildError error")]
    ThreadPoolBuildError(#[from] rayon::ThreadPoolBuildError),

    #[error("ProcessError")]
    ProcessError(#[from] paraseq::parallel::ProcessError),

    #[error("ParaseqFastaError")]
    ParaseqFastaError(#[from] paraseq::fasta::Error),

    #[error("ParaseqFastqError")]
    ParaseqFastqError(#[from] paraseq::fastq::Error),

    #[error(transparent)]
    Other(#[from] anyhow::Error),

    #[error("UTF8 error")]
    Utf(#[from] std::str::Utf8Error),

    #[error("Empty file: {0}")]
    EmptyFile(String),

    #[error("Invalid phred value")]
    InvalidPhredValue,

    #[error("Invalid figure types")]
    InvalidFigureType,
}
