mod bam2fastq;
mod dnbseqreorder;
mod fastqcheck;
mod illuminareorder;
mod namegroup;
mod readgroup;

#[derive(Debug, Clone, clap::Subcommand)]
pub enum Commands {
    Bam2fastq(bam2fastq::Fastq),
    FastqCheck(fastqcheck::FastqCheck),
    IlluminaFastqReorder(illuminareorder::IlluminaFastqReorder),
    DnbseqFastqReorder(dnbseqreorder::DNBSeqFastqReorder),
    PrintReadGroup(readgroup::CountReadGroup),
    NameGroup(namegroup::NameGroup),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::Bam2fastq(x) => x.run(),
            Commands::FastqCheck(x) => x.run(),
            Commands::IlluminaFastqReorder(x) => x.run(),
            Commands::DnbseqFastqReorder(x) => x.run(),
            Commands::PrintReadGroup(x) => x.run(),
            Commands::NameGroup(x) => x.run(),
        }
    }
}
