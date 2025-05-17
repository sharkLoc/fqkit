use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use anyhow::Error;
use log::{error, info};
use lowcharts::plot;
use plotters::prelude::*;
use paraseq::{
    fastq,
    fastx::Record,
    parallel::{ParallelProcessor, ParallelReader, ProcessError},
};
use parking_lot::Mutex;
use std::collections::HashMap;
use std::sync::Arc;

type Gctype = HashMap<u64, usize>;

#[derive(Clone)]
struct Gchash {
    thread_hash: Gctype,
    glob_hash: Arc<Mutex<Gctype>>,
}

impl Gchash {
    pub fn new() -> Self {
        Gchash {
            thread_hash: HashMap::new(),
            glob_hash: Arc::new(Mutex::new(HashMap::new())),
        }
    }
}

impl ParallelProcessor for Gchash {
    fn process_record<Rf: Record>(&mut self, record: Rf) -> Result<(), ProcessError> {
        let gc_count = record
            .seq()
            .iter()
            .filter(|x| *x == &b'G' || *x == &b'C')
            .count();
        let gc_ratio = (gc_count as f64 / record.seq().len() as f64 * 100.0).round() as u64;
        *self.thread_hash.entry(gc_ratio).or_insert(0) += 1usize;

        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<(), ProcessError> {
        let mut global_hash = self.glob_hash.lock();
        for (k, v) in self.thread_hash.iter() {
            global_hash.entry(*k).and_modify(|e| *e += *v).or_insert(*v);
        }

        // reset for next batch
        self.thread_hash = HashMap::new();
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub fn gc_content(
    fqin: Option<&String>,
    output: Option<&String>,
    show: bool,
    prefix: String,
    width: usize,
    height: usize,
    ylim: usize,
    types: &str,
    ncpu: usize,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), Error> {
    let fq_reader = file_reader(fqin).map(fastq::Reader::new)?;

    let gc_hash = Gchash::new();
    fq_reader.process_parallel(gc_hash.clone(), ncpu)?;
    let df_hash = gc_hash.glob_hash.lock();

    let mut fo = file_writer(output, compression_level, stdout_type)?;
    fo.write_all("GC(%)\tReads\tRatio(%)\n".as_bytes())?;
    let mut df_ret = vec![]; // data for PNG / SVG
    let mut df_num = vec![]; // data for histogram in terminal
    let total = df_hash.values().sum::<usize>() as f32;

    for i in 0..=100 {
        let num = *df_hash.get(&i).unwrap_or(&0);
        let v = (num as f32 * 10000.0 / total).round() / 100.0;
        df_ret.push(v);
        fo.write_all(format!("{}\t{}\t{}\n", i, num, v).as_bytes())?;
        if show {
            for _ in 0..num {
                df_num.push(i as f64)
            }
        }
    }
    fo.flush()?;

    //plot_gc(df_ret, prefix, width, height, ylim, types, quiet)?;
    if show {
        info!(
            "{}",
            plot::Histogram::new(
                &df_num,
                plot::HistogramOptions {
                    intervals: 20,
                    ..Default::default()
                }
            )
        );
    }
    plot_gc(df_ret, prefix, width, height, ylim, types)?;

    Ok(())
}

// GC content line plot
fn plot_gc(
    data: Vec<f32>,
    prefix: String,
    width: usize,
    height: usize,
    ylim: usize,
    types: &str,
) -> Result<(), Error> {
    if !["svg", "png"].contains(&types) {
        error!("{}", FqkitError::InvalidFigureType);
        std::process::exit(1);
    }
    if ylim > 100 {
        error!("invalid args ylim.");
        std::process::exit(1);
    }
    let name = if types == "png" {
        format!("{}.png", prefix)
    } else {
        format!("{}.svg", prefix)
    };
    info!("output gc content plot: {}", name);

    if types == "png" {
        let png = BitMapBackend::new(&name, (width as u32, height as u32)).into_drawing_area();
        png.fill(&WHITE)?;

        let mut charts = ChartBuilder::on(&png)
            .margin(10)
            .caption("GC distrbution plot", ("sans-serif", 30).into_font())
            .x_label_area_size(40)
            .y_label_area_size(40)
            //.build_cartesian_2d(0..100, 0..ylim)?;
            .build_cartesian_2d(0.0..100f32, 0.0..ylim as f32)?;

        charts
            .configure_mesh()
            .x_labels(20)
            .x_desc("GC CONTENT")
            .x_label_formatter(&|x| format!("{:.0}%", x))
            .y_labels(20)
            .y_label_formatter(&|x| format!("{:.0}%", x))
            .y_desc("percent")
            .draw()?;

        charts
            .draw_series(
                AreaSeries::new(
                    (0..).zip(data.iter()).map(|(x, y)| (x as f32, *y)),
                    0.,
                    GREEN.mix(0.5),
                )
                .border_style(GREEN),
            )?
            .label("GC content")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN));

        //charts.draw_series( LineSeries::new((0..).zip(data.iter()).map(|(x,y)| (x as f32, *y)), &BLACK) )?;

        charts
            .configure_series_labels()
            .background_style(WHITE.mix(0.9))
            .border_style(BLACK)
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;
    } else {
        let svg = SVGBackend::new(&name, (width as u32, height as u32)).into_drawing_area();
        svg.fill(&WHITE)?;

        let mut charts = ChartBuilder::on(&svg)
            .margin(10)
            .caption("GC distrbution plot", ("sans-serif", 30).into_font())
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.0..100f32, 0.0..ylim as f32)?;

        charts
            .configure_mesh()
            .x_labels(20)
            .x_desc("GC CONTENT")
            .x_label_formatter(&|x| format!("{:.0}%", x))
            .y_labels(20)
            .y_label_formatter(&|x| format!("{:.0}%", x))
            .y_desc("percent")
            .draw()?;

        charts
            .draw_series(
                AreaSeries::new(
                    (0..).zip(data.iter()).map(|(x, y)| (x as f32, *y)),
                    0.,
                    GREEN.mix(0.5),
                )
                .border_style(GREEN),
            )?
            .label("GC content")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN));

        charts
            .configure_series_labels()
            .background_style(WHITE.mix(0.9))
            .border_style(BLACK)
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;
    }
    Ok(())
}
