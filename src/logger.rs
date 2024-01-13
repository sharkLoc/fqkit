use std::{fs::File, io::{BufWriter, Write}};
use anyhow::{Ok, Result};
use chrono::Local;
use env_logger::{Builder,fmt::Color, Target};
use log::{LevelFilter,Level};


pub fn logger(
    verbose: String,
    logfile: &Option<&str>,
    quiet: bool,
) -> Result<(), anyhow::Error>{

    let mut level =  if verbose == "error".to_string() {
        LevelFilter::Error
    } else if verbose == "warn".to_string() {
        LevelFilter::Warn
    }else if verbose == "info".to_string() {
        LevelFilter::Info
    } else if verbose == "debug".to_string() {
        LevelFilter::Debug
    } else if verbose == "trace".to_string() {
        LevelFilter::Trace
    } else {
        LevelFilter::Off
    };
    if quiet {
        level = LevelFilter::Off;
    }

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
        writeln!(buf,
            "[{} {} - {}] {}",
            Local::now().format("%Y-%m-%dT%H:%M:%S"),
            style.value(record.level()),
            buf.style().set_color(Color::Rgb(90, 150, 150)).value(record.target()),
            record.args()
        )
    });
    
    // write log message in stderr or a file
    if let Some(file) = logfile {
        builder.target(Target::Pipe(log_writer(file)?)).filter(None, level).init();
    } else {
        builder.filter(None, level).init();
    }
    
    Ok(())
}


fn log_writer(file_out: &str) -> Result<Box<dyn Write + Send>> {
        let fp = File::create(file_out)?;
        Ok(Box::new(BufWriter::with_capacity(1, fp)))
}