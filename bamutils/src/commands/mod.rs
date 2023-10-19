mod concat_and_unify_read_group;
mod random_sampling;
mod rename_readname;
mod sequencing_error;

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Commands {
    SequencingError(sequencing_error::SequencingError),
    RandomSampling(random_sampling::RandomSampling),
    RenameReadname(rename_readname::RenameReadname),
    ConcatAndUnifyReadGroup(concat_and_unify_read_group::ConcatAndUnifyReadGroup),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::SequencingError(x) => x.run(),
            Commands::RandomSampling(x) => x.run(),
            Commands::RenameReadname(x) => x.run(),
            Commands::ConcatAndUnifyReadGroup(x) => x.run(),
        }
    }
}
