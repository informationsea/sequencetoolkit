mod random_sampling;
mod rename_readname;
mod sequencing_error;

use sequencetoolkit_common::Command;

pub(crate) const COMMANDS: &[&dyn Command] = &[
    &sequencing_error::SequencingError,
    &random_sampling::RandomSampling,
    &rename_readname::RenameReadname,
];
