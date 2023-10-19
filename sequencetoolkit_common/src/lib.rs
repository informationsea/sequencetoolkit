use clap::Args;

pub trait Command: Args {
    fn run(&self) -> anyhow::Result<()>;
}
