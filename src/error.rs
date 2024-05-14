use std::io;
use thiserror::Error;

#[derive(Debug,Error)]
pub enum FqkitError {
    #[error("stdin not detected")]
    StdinNotDetected,

    #[error("failed to open file: {0}")]
    IoError(#[from] io::Error),
    
    #[error("invalid output dir: {0}")]
    InvalidOutputDir(String),
    
    #[error("empty file: {0}")]
    EmptyFile(String),
    
    #[error("invalid phred value")]
    InvalidPhredValue,

}