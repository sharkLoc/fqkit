use anyhow::{Ok, Result};
use chrono::Local;
use env_logger::{Builder, Target, fmt::Color};
use log::{Level, LevelFilter};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

pub fn logger<P: AsRef<Path>>(
    verbose: u8,
    logfile: Option<P>,
    quiet: bool,
) -> Result<(), anyhow::Error> {
    let level = if quiet {
        LevelFilter::Off
    } else {
        match verbose {
            1 => LevelFilter::Error,
            2 => LevelFilter::Warn,
            3 => LevelFilter::Info,
            4 => LevelFilter::Debug,
            5 => LevelFilter::Trace,
            _ => LevelFilter::Off,
        }
    };

    let mut builder = Builder::from_default_env();
    builder.format(|buf, record| {
        let mut style = buf.style();
        match record.level() {
            Level::Error => {
                style.set_color(Color::Red).set_bold(true);
            }
            Level::Warn => {
                style.set_color(Color::Yellow).set_bold(true);
            }
            Level::Info => {
                style.set_color(Color::Green).set_bold(true);
            }
            Level::Debug => {
                style.set_color(Color::Blue).set_bold(true);
            }
            Level::Trace => {
                style.set_color(Color::Magenta).set_bold(true);
            }
        }

        writeln!(
            buf,
            "{} {} - {} {}",
            Local::now().format("%Y-%m-%dT%H:%M:%S"),
            style.value(record.level()),
            buf.style()
                .set_color(Color::Rgb(90, 150, 150))
                .value(record.target()),
            record.args()
        )
    });

    // write log message in stderr or a file
    if let Some(file) = logfile {
        builder
            .target(Target::Pipe(log_writer(file)?))
            .filter(None, level)
            .init();
    } else {
        builder.filter(None, level).init();
    }

    Ok(())
}

fn log_writer<P: AsRef<Path>>(file_out: P) -> Result<Box<dyn Write + Send>> {
    let fp = File::create(file_out)?;
    Ok(Box::new(BufWriter::with_capacity(1, fp)))
}
