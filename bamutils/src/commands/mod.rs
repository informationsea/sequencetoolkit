mod sequencing_error;

use sequencetoolkit_common::Command;

pub(crate) const COMMANDS: &[&dyn Command] = &[&sequencing_error::SequencingError];
