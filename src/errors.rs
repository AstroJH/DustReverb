use std::result;
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum Errors {
    #[error("Invalid IntervalResult")]
    FailedSearchInterval,
}

pub type Result<T> = result::Result<T, Errors>;