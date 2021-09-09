use std::convert::TryInto;
use std::error::Error;
use std::io::Read;
use std::io::Write;
use std::sync::{mpsc, Arc};

use clap::ArgMatches;
use crossbeam::queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressStyle};

use crate::config::{CHR_LENS, CHR_LENS_SMALL};
use crate::hmm;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let common_cells = hmm::get_cells(&sub_m)?;
    let num_common_cells = common_cells.len();
    info!(
        "Found {} cells in common assay, first cell={}",
        num_common_cells, common_cells[0]
    );

    let num_threads = 4;
    let num_chrs = CHR_LENS.len();
    let onlyone = sub_m.is_present("onlyone");
    let (num_states, chr_lens, chrs) = match onlyone {
        true => (2, CHR_LENS_SMALL, (0..1)),
        false => (12, CHR_LENS, (0..num_chrs)),
    };
    info!("Found total {} chromosomes", chr_lens.len());

    let in_path = carina::file::file_path_from_clap(&sub_m, "in_directory").unwrap();
    info!("Found input directory path: {:?}", in_path);

    let out_path = carina::file::file_path_from_clap(&sub_m, "out_directory").unwrap();
    info!("Found output directory path: {:?}", out_path);

    info!("Starting to read");
    //(0..num_chrs).rev().for_each(|chr_id| {
    chrs.rev().for_each(|chr_id| {
        let chr_name = format!("chr{}", chr_id+1);
        let num_bins = (chr_lens[chr_id] / 200) + 1;
        info!("Working on {}", chr_name);

        let pbar = ProgressBar::new(num_common_cells as u64);
        pbar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
                )
                .progress_chars("╢▌▌░╟"),
        );

        let q = Arc::new(ArrayQueue::<usize>::new(num_common_cells));
        (0..num_common_cells).for_each(|x| q.push(x).unwrap());
        let (tx, rx) = mpsc::sync_channel(num_threads);

        let chr_path = in_path.join(&chr_name);
        std::fs::create_dir_all(&chr_path).unwrap();

        let arc_in_path = Arc::new(&chr_path);
        let arc_common_cells = Arc::new(&common_cells);

        crossbeam::scope(|scope| {
            for _ in 0..num_threads {
                let tx = tx.clone();
                let reader = Arc::clone(&q);
                let arc_in_path = Arc::clone(&arc_in_path);
                let arc_common_cells = Arc::clone(&arc_common_cells);

                scope.spawn(move |_| loop {
                    match reader.pop() {
                        Some(cell_id) => {
                            let cell_file = arc_in_path.join(format!("{}.bin", arc_common_cells[cell_id]));
                            let mut file_handle = carina::file::bufreader_from_filepath(cell_file).unwrap();

                            let mut mat_u8_sizes = vec![0 as u8; 12];
                            file_handle.read_exact(&mut mat_u8_sizes).expect("can't read sizes");
                            let nnz: u32 = u32::from_le_bytes(mat_u8_sizes[0..4].try_into().unwrap());
                            let _nrows: u32 = u32::from_le_bytes(mat_u8_sizes[4..8].try_into().unwrap());
                            let _ncols: u32 = u32::from_le_bytes(mat_u8_sizes[8..12].try_into().unwrap());

                            let mut indices = vec![0 as u8; nnz as usize * 4];
                            let mut states = vec![0 as u8; nnz as usize];
                            let mut probs = vec![0 as u8; nnz as usize];

                            file_handle.read_exact(&mut probs).expect("can't read probs");
                            file_handle.read_exact(&mut states).expect("can't read states");
                            file_handle.read_exact(&mut indices).expect("can't read indices");

                            let states_size = states.len();
                            let mut state_indices = vec![Vec::<u8>::new(); num_states];
                            let mut state_probs = vec![Vec::new(); num_states];
                            for (idx, state) in states.into_iter().rev().enumerate() {
                                let idx = states_size - idx - 1;
                                let state: usize = u8::from_le(state) as usize;
                                state_probs[state].push(probs[idx]);
                                state_indices[state].extend(&indices[idx*4..(idx+1)*4]);
                            }

                            tx.send(Some((state_indices, state_probs, cell_id)))
                                .expect("Could not send mid data!");
                        }
                        None => {
                            tx.send(None).expect("Could not send end data!");
                            break;
                        }
                    }
                });
            }

            let chr_path = out_path.join(&chr_name);
            std::fs::create_dir_all(&chr_path).unwrap();

            let mut file_handles: Vec<std::io::BufWriter<std::fs::File>> = (0..num_states).map(|x| {
                let file_path = chr_path.join(&format!("{}.bin", x+1));
                std::io::BufWriter::new(std::fs::File::create(file_path).unwrap())
            }).collect();
            let mut cell_id_handle = std::io::BufWriter::new(std::fs::File::create(chr_path.join(&"cells.txt")).unwrap());

            let mut running_sums: Vec<u32> = vec![0; num_states];
            let mut sizes: Vec<Vec<u32>> = vec![vec![0]; num_states];
            let mut bin_probs: Vec<Vec<u8>> = vec![Vec::new(); num_states];
            let mut bin_indices: Vec<Vec<u8>> = vec![Vec::new(); num_states];
            let mut dead_thread_count = 0;
            for out_data in rx.iter() {
                match out_data {
                    Some((mut state_indices, mut state_probs, cell_id)) => {
                        for i in 0..num_states {
                            running_sums[i] += state_probs[i].len() as u32;
                            sizes[i].push(running_sums[i]);

                            bin_probs[i].append(&mut state_probs[i]);
                            bin_indices[i].append(&mut state_indices[i]);
                        }
                        writeln!(cell_id_handle, "{}", common_cells[cell_id]).unwrap();
                        pbar.inc(1);
                    } // end-Some
                    None => {
                        dead_thread_count += 1;
                        if dead_thread_count == num_threads {
                            drop(tx);
                            break;
                        }
                    } // end-None
                } // end-match
            } // end-for

            for i in 0..num_states {
                if bin_probs[i].len() == 0 { continue; }
                file_handles[i].write_all(&(num_bins as u32).to_le_bytes()).unwrap();
                file_handles[i].write_all(&(sizes[i].len() as u32).to_le_bytes()).unwrap();

                let bin_sizes: Vec<u8> = sizes[i].iter().map(|x| x.to_le_bytes()).collect::<Vec<[u8; 4]>>().concat();
                file_handles[i].write_all(&bin_sizes).unwrap();
                file_handles[i].write_all(&bin_probs[i]).unwrap();
                file_handles[i].write_all(&bin_indices[i]).unwrap();
            }
        })
        .unwrap(); //end crossbeam
        pbar.finish();
    }); // end for loop over chromosomes

    info!("All Done");
    Ok(())
}
