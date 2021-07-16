use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::ArgMatches;

use crossbeam::queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressStyle};
use std::sync::{mpsc, Arc};

pub fn get_moransi(
    weights: &sce::SingleCellExperiment<f32>,
    values: &sce::SingleCellExperiment<f32>,
    row_index: usize,
    row_sums: &[f32],
) -> Result<f32, Box<dyn Error>> {
    let n = values.cols();
    let val_it = values.counts().outer_iterator().nth(row_index).unwrap();

    let mut x_sum = 0.0;
    let mut x = vec![0.0_f32; n];
    for (col_ind, &val) in val_it.iter() {
        x_sum += val;
        x[col_ind] = val;
    }

    if x_sum == 0.0 {
        return Ok(0.0);
    }

    let x_mean = x_sum as f32 / n as f32;
    let z: Vec<f32> = x.iter().map(|&x| x as f32 - x_mean).collect();
    let v: f32 = z.iter().map(|x| x * x).sum();

    let mut w = 0.0;
    let mut cv = 0.0;
    for (i, row_iter) in weights.counts().outer_iterator().enumerate() {
        for (j, &wt) in row_iter.iter() {
            let norm_wt = wt / row_sums[i];
            w += norm_wt;
            cv += norm_wt * (z[i] as f32) * (z[j] as f32);
        }
    }

    Ok((n as f32 / w) * (cv / v))
}

pub fn get_gearyc(
    weights: &sce::SingleCellExperiment<f32>,
    values: &sce::SingleCellExperiment<f32>,
    row_index: usize,
    row_sums: &[f32],
) -> Result<f32, Box<dyn Error>> {
    let n = values.cols();
    let val_it = values.counts().outer_iterator().nth(row_index).unwrap();

    let mut x_sum = 0.0;
    let mut x = vec![0.0_f32; n];
    for (col_ind, &val) in val_it.iter() {
        x_sum += val;
        x[col_ind] = val;
    }

    if x_sum == 0.0 {
        return Ok(1.0);
    }

    let x_mean = x_sum as f32 / n as f32;
    let z: Vec<f32> = x.iter().map(|&x| x as f32 - x_mean).collect();
    let v: f32 = z.iter().map(|x| x * x).sum();

    let mut w = 0.0;
    let mut cv = 0.0;
    for (i, row_iter) in weights.counts().outer_iterator().enumerate() {
        for (j, &wt) in row_iter.iter() {
            let norm_wt = wt / row_sums[i];
            w += norm_wt;

            let x_mult = z[i] as f32 - z[j] as f32;
            cv += norm_wt * x_mult * x_mult;
        }
    }

    Ok(((n as f32 - 1.0) / (2.0 * w)) * (cv / v))
}

pub fn process(
    weights: &sce::SingleCellExperiment<f32>,
    values: &sce::SingleCellExperiment<f32>,
    mut ofile: BufWriter<File>,
    is_moransi: bool,
) -> Result<(), Box<dyn Error>> {
    let num_values = values.rows();
    let pbar = ProgressBar::new(num_values as u64);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    let row_sums: Vec<f32> = weights
        .counts()
        .outer_iterator()
        .map(|x| {
            let rsum: f32 = x.iter().map(|x| x.1).sum();
            if rsum == 0.0 {
                1.0
            } else {
                rsum
            }
        })
        .collect();

    let num_threads = 10;
    let q = Arc::new(ArrayQueue::<usize>::new(num_values));
    let arc_row_sums = Arc::new(row_sums);
    (0..num_values).for_each(|x| q.push(x).unwrap());

    let (tx, rx) = mpsc::sync_channel(num_threads);
    crossbeam::scope(|scope| {
        for _worker in 0..num_threads {
            let tx = tx.clone();
            let reader = Arc::clone(&q);
            let row_sums = Arc::clone(&arc_row_sums);

            scope.spawn(move |_| loop {
                match reader.pop() {
                    Some(index) => {
                        let stats = match is_moransi {
                            true => get_moransi(&weights, &values, index, &row_sums)
                                .expect("can't process rows"),
                            false => get_gearyc(&weights, &values, index, &row_sums)
                                .expect("can't process rows"),
                        };
                        tx.send(Some((index, stats)))
                            .expect("Could not send mid data!");
                    }
                    None => {
                        tx.send(None).expect("Could not send end data!");
                        break;
                    }
                }
            });
        }

        let mut dead_thread_count = 0;
        let column_names = values.col_names();
        for out_data in rx.iter() {
            match out_data {
                Some((index, stats)) => {
                    pbar.inc(1);
                    writeln!(ofile, "{}\t{}", column_names[index], stats)
                        .expect("can't write to file");
                } // end-Some
                None => {
                    dead_thread_count += 1;
                    if dead_thread_count == num_threads {
                        drop(tx);

                        // consume what's remaining
                        for out_data in rx.iter() {
                            pbar.inc(1);
                            out_data.map_or((), |(index, stats)| {
                                writeln!(ofile, "{}\t{}", column_names[index], stats)
                                    .expect("can't write to file");
                            });
                        }

                        break;
                    } // end if
                } // end-None
            } // end-match
        } // end-for
    })
    .unwrap(); //end crossbeam

    pbar.finish();
    Ok(())
}

pub fn generate_stats(
    weights_file_path: PathBuf,
    values_file_path: PathBuf,
    ofile: BufWriter<File>,
    method: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    let wt_mat: sce::SingleCellExperiment<f32> =
        sce::SingleCellExperiment::from_tenx_v2(weights_file_path)?;
    println!("Weights: {:?}", wt_mat);

    let val_mat: sce::SingleCellExperiment<f32> =
        sce::SingleCellExperiment::from_tenx_v2(values_file_path)?;
    println!("Values: {:?}", val_mat);

    match method {
        Some("Moransi") => {
            info!("Starting Moran's I");
            process(&wt_mat, &val_mat, ofile, true)?;
        }
        Some("Gearyc") => {
            info!("Starting Moran's I");
            process(&wt_mat, &val_mat, ofile, false)?;
        }
        _ => unreachable!(),
    };

    Ok(())
}

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let weights_file_path = carina::file::file_path_from_clap(sub_m, "weights")?;
    let values_file_path = carina::file::file_path_from_clap(sub_m, "values")?;
    let ofile = carina::file::bufwriter_from_clap(sub_m, "output")?;

    generate_stats(
        weights_file_path,
        values_file_path,
        ofile,
        sub_m.value_of("method"),
    )?;

    info!("All done");
    Ok(())
}
