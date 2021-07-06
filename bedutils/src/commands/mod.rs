mod extend_bed;
mod merge_bed;

use crate::Command;

pub(crate) const COMMANDS: &[&dyn Command] = &[&merge_bed::MergeBed, &extend_bed::ExtendBed];
