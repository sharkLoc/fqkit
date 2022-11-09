use std::io::BufRead;
use std::collections::BTreeMap;
use plotters::prelude::*;
use colored::*;

use crate::utils::file_reader;


// get cycle result
pub fn cycle_data(file: &Option<&str>) -> std::io::Result<Vec<BTreeMap<usize,f64>>>{
    let mut cyc: Vec<BTreeMap<usize,f64>> = Vec::new();
    let fp = file_reader(file)?;
    for _ in 0..5 {  cyc.push(BTreeMap::new()); }

    for (idx, line) in fp.lines().flatten().enumerate() {
        if idx == 0 { continue; }
        let mut info = line.split_whitespace().take(6).collect::<Vec<&str>>();
        info.remove(0);
        for (i,e) in info.iter().enumerate() {
            let rate = e.strip_suffix("%)").unwrap().split("(").collect::<Vec<&str>>()[1].parse::<f64>().expect("error: not a number");
            cyc[i].insert(idx, rate);
        }
    }
    Ok(cyc)
}

// line plot for base A T G C N rate in position
pub fn plot_line(
    data: Vec<BTreeMap<usize,f64>>, 
    prefix: String, 
    width: usize, 
    height: usize,
    ylim: f32,
    types: &str) -> Result<(), Box<dyn std::error::Error>> {
    if !["svg", "png"].contains(&types) {
        eprintln!("{}","[error]: invalid args types.".red());
        std::process::exit(1);
    }
    let max_len = *data[0].iter().last().unwrap().0 as f32;
    if ylim < 0.0 {
        eprintln!("{}","[error]: invalid args ylim.".red());
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
    Ok(())
}