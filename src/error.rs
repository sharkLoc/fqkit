use std::io;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FqkitError {
    #[error("Stdin not detected")]
    StdinNotDetected,

    #[error("Failed to open file: {0}")]
    IoError(#[from] io::Error),

    #[error("Invalid output dir: {0}")]
    InvalidOutputDir(String),

    #[error("ThreadPoolBuildError error")]
    ThreadPoolBuildError(#[from] rayon::ThreadPoolBuildError),

    #[error("Empty file: {0}")]
    EmptyFile(String),

    #[error("Invalid phred value")]
    InvalidPhredValue,

    #[error("Invalid figure types")]
    InvalidFigureType,
}
