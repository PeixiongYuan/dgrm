
use thiserror::Error;

#[derive(Error, Debug)]
pub enum DgrmError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("CSV parsing error: {0}")]
    Csv(#[from] csv::Error),

    #[error("Invalid data format: {0}")]
    InvalidFormat(String),

    #[error("Matrix dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },

    #[error("File not found: {path}")]
    FileNotFound { path: String },
}

pub type Result<T> = std::result::Result<T, DgrmError>;