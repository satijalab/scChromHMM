use crate::config::ProbT;
use crate::config::{CHR_LENS, CHR_LENS_SMALL, WINDOW_SIZE};
use crate::fragment::Fragment;
use crate::model;
use crate::quantify;
use crate::record::{AssayRecords, Experiment};
use bio::data_structures::interval_tree::IntervalTree;

use clap::ArgMatches;
use crossbeam::queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressStyle};
use std::sync::{mpsc, Arc};

use std::collections::HashMap;
use std::error::Error;
use std::io::BufRead;
use std::io::Write;
use std::ops::Range;

//use flate2::write::GzEncoder;
//use flate2::Compression;

pub fn get_cells(sub_m: &ArgMatches) -> Result<Vec<String>, Box<dyn Error>> {
    // reading in cell names of common assay.
    let cells_file_path = carina::file::file_path_from_clap(sub_m, "common_cells")?;
    let cells_reader = carina::file::bufreader_from_filepath(cells_file_path)?;

    let mut common_cells: Vec<String> = Vec::with_capacity(5000);
    for line in cells_reader.lines() {
        common_cells.push(line.unwrap());
    }

    Ok(common_cells)
}

fn get_anchors(
    sub_m: &ArgMatches,
    common_cells: &[String],
) -> Result<Vec<HashMap<u64, HashMap<u32, ProbT>>>, Box<dyn Error>> {
    // reading anchors
    let string_index_common_cells: HashMap<String, u32> = common_cells
        .iter()
        .enumerate()
        .map(|(i, x)| (x.clone(), i as u32))
        .collect();

    let mut vec_anchor_triplets = Vec::with_capacity(5);

    let anchor_file_paths = carina::file::files_path_from_clap(sub_m, "anchors")?;
    for file_path in anchor_file_paths {
        let mut anchor_triplets: HashMap<u64, HashMap<u32, ProbT>> = HashMap::with_capacity(10_000);

        let reader = carina::file::bufreader_from_filepath(file_path)?;
        for line in reader.lines() {
            let record = line.unwrap();
            let toks: Vec<&str> = record.split_whitespace().collect();

            let common_cell_index: u32 = *string_index_common_cells
                .get(toks[0])
                .expect("can't find cell name in the common cell list");

            let cb_str: Vec<&str> = toks[1].split('-').collect();
            let cb_id = cb_str.get(1).unwrap().parse::<u8>().unwrap();

            let assay_cell_index: u64 =
                carina::barcode::cb_string_to_u64_with_id(cb_str[0].as_bytes(), 16, cb_id)?;

            let score = toks[2].parse::<ProbT>().unwrap();
            anchor_triplets
                .entry(assay_cell_index)
                .or_insert_with(HashMap::new)
                .insert(common_cell_index, score);
        }

        vec_anchor_triplets.push(anchor_triplets);
    }

    Ok(vec_anchor_triplets)
}

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let num_threads: usize = sub_m.value_of("threads").unwrap().parse().unwrap();

    let hmm = model::get_hmm_params(&sub_m)?;
    info!("Read HMM model paramers: {:?}", hmm);

    let common_cells = get_cells(&sub_m)?;
    let num_common_cells = common_cells.len();
    info!(
        "Found {} cells in common assay, first cell={}",
        common_cells.len(),
        common_cells[0]
    );

    let vec_anchor_triplets = get_anchors(&sub_m, &common_cells)?;
    let num_assays = vec_anchor_triplets.len();
    let assay_num_cells: Vec<usize> = vec_anchor_triplets.iter().map(|x| x.len()).collect();
    let assay_num_anchors: Vec<usize> = vec_anchor_triplets
        .iter()
        .map(|x| x.iter().map(|(_, z)| z.len()).sum())
        .collect();
    info!(
        "Found {} assays with {:?} cells and {:?} anchors",
        num_assays, assay_num_cells, assay_num_anchors
    );

    let fragment_file_paths = carina::file::files_path_from_clap(sub_m, "fragments")?;
    let mut frags: Vec<Fragment> = fragment_file_paths
        .into_iter()
        .map(Fragment::from_pathbuf)
        .collect();

    let num_chrs = CHR_LENS.len();
    let onlyone = sub_m.is_present("onlyone");

    let (chrs, chr_lens) = match onlyone {
        true => {
            ((0..1), CHR_LENS_SMALL)
        },
        false => {
            ((0..num_chrs), CHR_LENS)
        },
    };
    info!("Found total {} chromosomes", chrs.len());

    info!("Starting forward backward");
    //(0..num_chrs).rev().take(1).for_each(|chr_id| {
    chrs.rev().for_each(|chr_id| {
        let chr_name = format!("chr{}", chr_id+1);
        let tids: Vec<u64> = frags.iter().map(|x| x.tid(&chr_name)).collect();
        info!("Working on {}", chr_name);

        let range = Range {
            start: 0,
            end: chr_lens[chr_id],
        };
        let assay_data: Vec<AssayRecords<ProbT>> = frags
            .iter_mut()
            .enumerate()
            .map(|(i, x)| {
                let cell_records = x.fetch(
                    tids[i],
                    &range,
                    &vec_anchor_triplets.get(i).unwrap(),
                    num_common_cells,
                );

                AssayRecords::new(cell_records)
            })
            .collect();
        let exp = Experiment::new(assay_data);

        let pbar = ProgressBar::new(num_common_cells as u64);
        pbar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
                )
                .progress_chars("╢▌▌░╟"),
        );

        let only_imputed = false;
        if only_imputed {
            let mut files: Vec<std::io::BufWriter<std::fs::File>> = (0..num_assays).map(|id| {
                let path = std::path::Path::new("/mnt/scratch1/avi/Indus/data/out").join(&chr_name).join(format!("{}.txt", id));
                let f = std::fs::File::create(path).unwrap();
                std::io::BufWriter::new(f)
            }).collect();

            (0..num_common_cells).for_each(|cell_id| {
                pbar.inc(1);
                let cell_data = exp.get_cell_data(cell_id);
                let itrees: Vec<IntervalTree<u32, ProbT>> = cell_data
                    .into_iter()
                    .map(|cell_records| {
                        let mut tree = IntervalTree::new();
                        for record in cell_records.records() {
                            tree.insert(record.range(), record.id());
                        }
                        tree
                    })
                    .collect();

                let get_obv_list = |start: usize, end: usize| {
                    let observation_list: Vec<Vec<ProbT>> = (start..end)
                        .step_by(WINDOW_SIZE)
                        .map(|qstart| {
                            let qrange = qstart as u32..(qstart + WINDOW_SIZE) as u32;
                            let cts: Vec<ProbT> = itrees
                                .iter()
                                .map(|tree| {
                                    let vals: Vec<ProbT> = tree.find(&qrange).map(|x| *x.data()).collect();
                                    vals.iter().sum()
                                })
                                .collect();
                            cts
                        })
                        .collect();
                    observation_list
                };

                let observation_list = get_obv_list(0, chr_lens[chr_id] as usize);
                for observation in observation_list {
                    (0..observation.len()).for_each(|id| {
                        if observation[id] != 0.0 {
                            writeln!(files[id], "{:.8}", observation[id]).unwrap();
                        }
                    });
                }
            });
        } else {
            let out_path =
                std::path::Path::new(sub_m.value_of("output").unwrap()).join(&chr_name);
            std::fs::create_dir_all(&out_path).unwrap();

            let q = Arc::new(ArrayQueue::<usize>::new(num_common_cells));
            //(0..num_common_cells).filter(|&x| x == 2840).for_each(|x| q.push(x).unwrap());
            (0..num_common_cells).for_each(|x| q.push(x).unwrap());

            let (tx, rx) = mpsc::sync_channel(num_threads);

            let arc_hmm = Arc::new(&hmm);
            let arc_exp = Arc::new(&exp);
            let arc_out_path = Arc::new(&out_path);
            let arc_common_cells = Arc::new(&common_cells);

            let num_states = hmm.num_states();
            let chr_len = chr_lens[chr_id] as usize;
            crossbeam::scope(|scope| {
                for _ in 0..num_threads {
                    let tx = tx.clone();
                    let reader = Arc::clone(&q);
                    let arc_hmm = Arc::clone(&arc_hmm);
                    let arc_exp = Arc::clone(&arc_exp);
                    let arc_out_path = Arc::clone(&arc_out_path);
                    let arc_common_cells = Arc::clone(&arc_common_cells);

                    let mut posterior = Vec::with_capacity(chr_len * num_states / 400);
                    let mut fprob = vec![vec![0.0; arc_hmm.num_states()]; (chr_len / 200) + 1];

                    scope.spawn(move |_| loop {
                        match reader.pop() {
                            Some(cell_id) => {
                                posterior.clear();
                                let cell_data = arc_exp.get_cell_data(cell_id);
                                quantify::run_fwd_bkw(cell_data, &arc_hmm, &mut fprob, &mut posterior, chr_len).unwrap();

                                let out_file = arc_out_path.join(format!("{}.bin", arc_common_cells[cell_id]));
                                let num_posteriors = posterior.len();
                                let mut bin_mat: Vec<u8> = vec![(
                                    num_posteriors as u32).to_le_bytes(),
                                    ((chr_len / 200 + 1) as u32).to_le_bytes(),
                                    (num_states as u32).to_le_bytes()
                                ].concat();

                                let mut indices = posterior.iter().map(|x| (x.0 as u32).to_le_bytes()).collect::<Vec<[u8; 4]>>().concat();
                                let mut state = posterior.iter().map(|x| (x.1 as u8).to_le_bytes()).collect::<Vec<[u8; 1]>>().concat();
                                let mut value = posterior.iter().map(|x| ( (x.2 * 100.0).round() as u8).to_le_bytes()).collect::<Vec<[u8; 1]>>().concat();
                                bin_mat.append(&mut value);
                                bin_mat.append(&mut state);
                                bin_mat.append(&mut indices);

                                tx.send(Some((bin_mat, out_file)))
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
                for out_data in rx.iter() {
                    match out_data {
                        Some((mat, out_file)) => {
                            pbar.inc(1);
                            write_binary(out_file, mat).unwrap();
                            //write_matrix_market(out_file, &mat.to_csr()).unwrap();
                        } // end-Some
                        None => {
                            dead_thread_count += 1;
                            if dead_thread_count == num_threads {
                                drop(tx);

                                for out_data in rx.iter() {
                                    pbar.inc(1);
                                    out_data.map_or((), |(mat, out_file)| {
                                        write_binary(out_file, mat).unwrap();
                                        //write_matrix_market(out_file, &mat.to_csr()).unwrap();
                                    });
                                }
                                break;
                            }
                        } // end-None
                    } // end-match
                } // end-for
            })
            .unwrap(); //end crossbeam
        }
        pbar.finish();
    });
    info!("All Done");

    Ok(())
}

pub fn write_binary(path: std::path::PathBuf, mat: Vec<u8>) -> Result<(), Box<dyn Error>> {
    let f = std::fs::File::create(path)?;
    //let mut file = GzEncoder::new(f, Compression::default());
    let mut file = std::io::BufWriter::new(f);
    //let mut file = std::io::BufWriter::new(snap::write::FrameEncoder::new(f));

    // entries
    file.write_all(&mat)?;

    Ok(())
}
