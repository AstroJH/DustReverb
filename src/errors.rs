use std::result;
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum Errors {
    #[error("Invalid IntervalResult.")]
    FailedSearchInterval,

    #[error("Error opening file.")]
    FailedOpenFile,

    #[error("The header of CSV file has error.")]
    CsvHeaderError,
}

pub type Result<T> = result::Result<T, Errors>;