mod create_db;
mod genome_position;
mod transcript_position;
use crate::Command;

pub(crate) const COMMANDS: &[&dyn Command] = &[
    &create_db::CreateDb,
    &genome_position::GenomePosition,
    &transcript_position::TranscriptPosition,
];
