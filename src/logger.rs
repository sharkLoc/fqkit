use anyhow::Ok;
use chrono::Local;
use env_logger::{Builder,fmt::Color};
use log::{LevelFilter,Level};
use std::io::Write;


pub fn logger(verbose: String) -> Result<(), anyhow::Error>{

    let level =  if verbose == "error".to_string() {
        LevelFilter::Error
    } else if verbose == "warn".to_string() {
        LevelFilter::Warn
    }else if verbose == "info".to_string() {
        LevelFilter::Info
    } else if verbose == "debug".to_string() {
        LevelFilter::Debug
    } else {
        LevelFilter::Trace
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
        writeln!(buf,
            "[{} {} - {}] {}",
            Local::now().format("%Y-%m-%dT%H:%M:%S"),
            style.value(record.level()),
            buf.style().set_color(Color::Rgb(90, 150, 150)).value(record.target()),
            record.args()
        )
    })
    .filter(None, level)
    .init();

    Ok(())
}