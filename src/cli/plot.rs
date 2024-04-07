use std::io::BufRead;
use std::collections::BTreeMap;
use plotters::prelude::*;
use log::*;
use textplots::{Chart, ColorPlot, LabelBuilder, LabelFormat, Shape};
use colored::*;
use std::time::Instant;
use crate::utils::file_reader;
use anyhow::Result;

// get cycle result
pub fn cycle_data(file: Option<&String>) -> Result<Vec<BTreeMap<usize,f64>>>{
    let mut cyc: Vec<BTreeMap<usize,f64>> = Vec::new();
    let fp = file_reader(file)?;
    for _ in 0..5 {  cyc.push(BTreeMap::new()); }

    for (idx, line) in fp.lines().map_while(std::io::Result::ok).enumerate() {
        if idx == 0 { continue; }
        let mut info = line.split_whitespace().take(6).collect::<Vec<&str>>();
        info.remove(0);
        for (i,e) in info.iter().enumerate() {
            let rate = e.strip_suffix("%)").unwrap().split('(').collect::<Vec<&str>>()[1].parse::<f64>().expect("error: not a number");
            cyc[i].insert(idx, rate);
        }
    }
    Ok(cyc)
}

// line plot for base A T G C N rate in position
pub fn plot_line(
    data: Vec<BTreeMap<usize,f64>>, 
    show: bool,
    prefix: String, 
    width: usize, 
    height: usize,
    ylim: f32,
    types: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();

    if !["svg", "png"].contains(&types) {
        error!("invalid args types.");
        std::process::exit(1);
    }
    let max_len = *data[0].iter().last().unwrap().0 as f32;
    if ylim < 0.0 {
        error!("invalid args ylim.");
        std::process::exit(1);
    }
    let name = if types == "png" {format!("{}.png",prefix)} else {format!("{}.svg",prefix)};
    if types == "png" {
        let png = BitMapBackend::new(&name, (width as u32, height as u32)).into_drawing_area();
        png.fill(&WHITE)?;

        let mut charts = ChartBuilder::on(&png)
            .margin(10)
            .caption("Base distrbution plot", ("sans-serif", 40).into_font())
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.1..max_len , -0.5f32..ylim)?;

        charts
            .configure_mesh()
            .x_labels(20)
            .x_desc("position")
            .x_label_formatter(&|x| format!("{:.0}",x))
            .y_labels(10)
            .y_label_formatter(&|x| format!("{:.1}", x))
            .y_desc("percent")
            .draw()?;

        let nt_a= data[0].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();
        let nt_t= data[1].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();
        let nt_g= data[2].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();
        let nt_c= data[3].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();
        let nt_n= data[4].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();

        if show {
            let reads_len = nt_a.len();
            println!( "read length: {}\t\t{}\t{}\t{}\t{}\t{}", 
                    reads_len, 
                    "A".truecolor(255, 0, 0).bold(),
                    "T".truecolor(0, 255, 0).bold(),
                    "G".truecolor(255, 255, 0).bold(),
                    "C".truecolor(0, 0, 255).bold(),
                    "N".truecolor(0, 255, 255).bold());
            Chart::new_with_y_range(200, 80, 0.0, reads_len as f32, 0.0, ylim)
                .linecolorplot(&Shape::Lines(nt_c.as_slice()), rgb::RGB { r: 0, g: 0, b: 255 }) // C blue
                .linecolorplot(&Shape::Lines(nt_t.as_slice()), rgb::RGB { r: 0, g: 255, b: 0 }) // T green
                .linecolorplot(&Shape::Lines(nt_a.as_slice()), rgb::RGB { r: 255, g: 0, b: 0 }) // A red
                .linecolorplot(&Shape::Lines(nt_g.as_slice()), rgb::RGB { r: 255, g: 255, b: 0 }) // G yellow
                .linecolorplot(&Shape::Lines(nt_n.as_slice()), rgb::RGB { r: 0, g: 255, b: 255 }) // N white
                //.x_label_format(LabelFormat::Custom(Box::new( |_| { String::from("position")})))
                .y_label_format(LabelFormat::Value)
                .nice();
        }

        charts
            .draw_series(LineSeries::new(nt_a, RED))
            .unwrap()
            .label("A")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], RED));

        charts
            .draw_series(LineSeries::new(nt_t, GREEN))
            .unwrap()
            .label("T")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], GREEN));
        
        charts
            .draw_series(LineSeries::new(nt_g, YELLOW))
            .unwrap()
            .label("G")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], YELLOW));
        
        charts
            .draw_series(LineSeries::new(nt_c, BLACK))
            .unwrap()
            .label("C")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], BLACK));
        
        charts
            .draw_series(LineSeries::new(nt_n, BLUE))
            .unwrap()
            .label("N")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], BLUE));

    
    
        charts
            .configure_series_labels()
            .background_style(&WHITE.mix(0.9))
            .border_style(&BLACK)
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;

    } else {
        let svg = SVGBackend::new(&name, (width as u32, height as u32)).into_drawing_area();
        svg.fill(&WHITE)?;

        let mut charts = ChartBuilder::on(&svg)
            .margin(10)
            .caption("Base distrbution plot", ("sans-serif", 40).into_font())
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.1..max_len , -0.5f32..ylim)?;

        charts
            .configure_mesh()
            .x_labels(20)
            .x_desc("position")
            .x_label_formatter(&|x| format!("{:.0}",x))
            .y_labels(10)
            .y_label_formatter(&|x| format!("{:.1}", x))
            .y_desc("percent")
            .draw()?;

        let nt_a= data[0].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();
        let nt_t= data[1].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();
        let nt_g= data[2].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();
        let nt_c= data[3].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();
        let nt_n= data[4].iter().map(|(k,v)|(*k as f32, *v as f32)).collect::<Vec<(f32,f32)>>();

        if show {
            let reads_len = nt_a.len();
            println!( "read length: {}\t\t{}\t{}\t{}\t{}\t{}", 
                    reads_len, 
                    "A".truecolor(255, 0, 0).bold(),
                    "T".truecolor(0, 255, 0).bold(),
                    "G".truecolor(255, 255, 0).bold(),
                    "C".truecolor(0, 0, 255).bold(),
                    "N".truecolor(0, 255, 255).bold());
            Chart::new_with_y_range(200, 80, 0.0, reads_len as f32, 0.0, ylim)
                .linecolorplot(&Shape::Lines(nt_c.as_slice()), rgb::RGB { r: 0, g: 0, b: 255 }) // C blue
                .linecolorplot(&Shape::Lines(nt_t.as_slice()), rgb::RGB { r: 0, g: 255, b: 0 }) // T green
                .linecolorplot(&Shape::Lines(nt_a.as_slice()), rgb::RGB { r: 255, g: 0, b: 0 }) // A red
                .linecolorplot(&Shape::Lines(nt_g.as_slice()), rgb::RGB { r: 255, g: 255, b: 0 }) // G yellow
                .linecolorplot(&Shape::Lines(nt_n.as_slice()), rgb::RGB { r: 0, g: 255, b: 255 }) // N white
                .y_label_format(LabelFormat::Value)
                .nice();
        }

        charts
            .draw_series(LineSeries::new(nt_a, RED))
            .unwrap()
            .label("A")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], RED));

        charts
            .draw_series(LineSeries::new(nt_t, GREEN))
            .unwrap()
            .label("T")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], GREEN));
        
        charts
            .draw_series(LineSeries::new(nt_g, YELLOW))
            .unwrap()
            .label("G")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], YELLOW));
        
        charts
            .draw_series(LineSeries::new(nt_c, BLACK))
            .unwrap()
            .label("C")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], BLACK));
        
        charts
            .draw_series(LineSeries::new(nt_n, BLUE))
            .unwrap()
            .label("N")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], BLUE));

        charts
            .configure_series_labels()
            .background_style(&WHITE.mix(0.9))
            .border_style(&BLACK)
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;

    } 
    
    info!("time elapsed is: {:?}",start.elapsed());
    Ok(())
}
