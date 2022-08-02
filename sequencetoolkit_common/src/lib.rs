use clap::{crate_authors, crate_version, App, ArgMatches, SubCommand};

pub trait Command {
    fn cli(&self) -> App<'static, 'static> {
        self.config_subcommand(SubCommand::with_name(self.command_name()))
            .version(crate_version!())
            .author(crate_authors!())
    }
    fn command_name(&self) -> &'static str;
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static>;
    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()>;
}
