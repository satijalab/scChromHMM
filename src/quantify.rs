use bio::data_structures::interval_tree::IntervalTree;

use crate::config::{ProbT, MIN_PROB, WINDOW_SIZE};
use crate::model::Hmm;
use crate::record::CellRecords;

use std::error::Error;

fn forward(
    observations: &[Vec<ProbT>],
    hmm: &Hmm,
    fprob: &mut Vec<Vec<ProbT>>,
    num_states: usize,
    num_observations: usize,
) -> ProbT {
    let mut f_prev: Vec<ProbT> = (0..num_states)
        .map(|state| hmm.get_emission_prob(state, &observations[0]) * hmm.get_init_prob(state))
        .collect();
    let prob_norm: ProbT = f_prev.iter().sum();
    f_prev.iter_mut().for_each(|x| *x /= prob_norm);
    fprob[0].clone_from(&f_prev);

    let mut f_curr = vec![0.0; num_states];
    for i in 1..num_observations {
        for (state, f_item) in f_curr.iter_mut().enumerate().take(num_states) {
            let mut prev_f_sum = 0.0;
            (0..num_states).for_each(|prev_state| {
                prev_f_sum += f_prev[prev_state] * hmm.get_transition_prob(prev_state, state);
            });
            *f_item = hmm.get_emission_prob(state, &observations[i]) * prev_f_sum;
        }

        let prob_norm: ProbT = f_curr.iter().sum();
        f_curr.iter_mut().for_each(|x| *x /= prob_norm);

        fprob[i].clone_from(&f_curr);
        f_prev.clone_from(&f_curr);
    }

    let mut norm = f_curr.into_iter().map(|x| x * 0.1).sum();
    if norm == 0.0 {
        norm = 1.0;
    }

    norm
}

#[inline]
fn update_triplet(
    i: usize,
    posterior: &mut Vec<(usize, usize, ProbT)>,
    b_curr: &[ProbT],
    num_states: usize,
    norm: ProbT,
    fprob: &[Vec<ProbT>],
) {
    let is_valid_state = |state: usize| match state {
        //0 | 1 | 2 | 4 | 10 | 11 => true,
        0 | 1 | 2 | 3 | 8 | 9 | 11 => true,
        //0 | 1  => true,
        _ => false,
    };

    let probs: Vec<ProbT> = (0..num_states)
        .map(|state| fprob[i][state] * b_curr[state] / norm)
        .collect();
    let state_norm: ProbT = probs.iter().sum();
    probs.into_iter().enumerate().for_each(|(state, prob)| {
        let prob = prob / state_norm;
        if (prob > MIN_PROB) & is_valid_state(state) {
            posterior.push((i, state, prob))
        }
    });
}

fn backward(
    observations: Vec<Vec<ProbT>>,
    hmm: &Hmm,
    norm: ProbT,
    num_states: usize,
    num_observations: usize,
    fprob: &[Vec<ProbT>],
    posterior: &mut Vec<(usize, usize, ProbT)>,
) {
    let mut b_curr = vec![0.1; num_states];
    let mut b_prev = vec![0.1; num_states];

    update_triplet(
        num_observations - 1,
        posterior,
        &b_curr,
        num_states,
        norm,
        fprob,
    );
    for i in (1..num_observations).rev() {
        let obv_emissions: Vec<ProbT> = (0..num_states)
            .into_iter()
            .map(|state| hmm.get_emission_prob(state, &observations[i]))
            .collect();

        b_curr.iter_mut().for_each(|i| *i = 0.0);
        for (state, b_item) in b_curr.iter_mut().enumerate().take(num_states) {
            for next_state in 0..num_states {
                *b_item += hmm.get_transition_prob(state, next_state)
                    * obv_emissions[next_state]
                    * b_prev[next_state];
            }
        }

        let prob_norm: ProbT = b_curr.iter().sum();
        b_curr.iter_mut().for_each(|x| *x /= prob_norm);

        b_prev.clone_from(&b_curr);
        update_triplet(i - 1, posterior, &b_curr, num_states, norm, fprob);
    }
}

fn get_posterior(
    observations: Vec<Vec<ProbT>>,
    hmm: &Hmm,
    fprob: &mut Vec<Vec<ProbT>>,
    posterior: &mut Vec<(usize, usize, ProbT)>,
) {
    let num_states = hmm.num_states();
    let num_assays = hmm.num_assays();
    let num_observations = observations.len();

    assert!(num_assays == observations[0].len());
    let norm = forward(&observations, &hmm, fprob, num_states, num_observations);

    backward(
        observations,
        &hmm,
        norm,
        num_states,
        num_observations,
        fprob,
        posterior,
    );
}

pub fn run_fwd_bkw(
    cell_records: Vec<&CellRecords<ProbT>>,
    hmm: &Hmm,
    fprob: &mut Vec<Vec<ProbT>>,
    posterior: &mut Vec<(usize, usize, ProbT)>,
    chr_len: usize,
) -> Result<(), Box<dyn Error>> {
    let itrees: Vec<IntervalTree<u32, ProbT>> = cell_records
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
                let qstart: usize = std::cmp::max(0, qstart as i32 - 1) as usize;
                let qrange = qstart as u32..(qstart + WINDOW_SIZE + 1) as u32;
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

    //let qstart = 132492 * 200;
    //let qrange = qstart as u32 - 1..(qstart + WINDOW_SIZE) as u32;
    //let cts: Vec<Vec<(u32, u32, f32)>> = itrees
    //    .iter()
    //    .map(|tree| {
    //        let vals: Vec<(u32, u32, f32)> = tree.find(&qrange).map(|x| (x.interval().start, x.interval().end, *x.data())).collect();
    //        vals
    //    })
    //    .collect();
    //println!("{:?}", qrange);
    //println!("{:?}", cts);

    let observation_list = get_obv_list(0, chr_len);
    //println!("{:?}", observation_list[132491]);
    //println!("{:?}", observation_list[132492]);
    //println!("{:?}", observation_list[132493]);
    //println!("{:?}", observation_list[132494]);

    get_posterior(observation_list, hmm, fprob, posterior);

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_fwd_bkw() {
        let path = std::path::PathBuf::from("/mnt/scratch1/avi/Indus/data/model_test.txt");
        let file_reader = carina::file::bufreader_from_filepath(path).unwrap();
        let hmm = crate::model::Hmm::new(file_reader);
        let mut fprob = vec![vec![0.0; 2]; 3];
        let mut post = Vec::new();

        crate::quantify::get_posterior(
            vec![
                vec![1.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
                vec![0.0, 0.0, 1.0],
            ],
            &hmm,
            &mut fprob,
            &mut post,
        );
        let mut mat = sprs::TriMat::new((3, 2));
        post.into_iter()
            .for_each(|(x, y, z)| mat.add_triplet(x, y, z));
        let probs: Vec<String> = mat
            .to_csr::<usize>()
            .data()
            .into_iter()
            .map(|&x| format!("{:.4}", x))
            .collect();

        assert_eq!(
            probs,
            ["0.9342", "0.0658", "0.6665", "0.3335", "0.1209", "0.8791"]
        );
    }

    #[test]
    fn test_fwd_bkw_full() {
        let path = std::path::PathBuf::from("/home/srivastavaa/parazodiac/Indus/data/model_12.txt");
        let file_reader = carina::file::bufreader_from_filepath(path).unwrap();
        let hmm = crate::model::Hmm::new(file_reader);
        let mut fprob = vec![vec![0.0; 12]; 5];
        let mut post = Vec::new();

        crate::quantify::get_posterior(
            vec![
                vec![1.0, 0.0, 0.0, 0.0, 0.0],
                vec![1.0, 0.0, 0.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0, 0.0, 0.0],
                vec![1.0, 0.0, 0.0, 0.0, 0.0],
            ],
            &hmm,
            &mut fprob,
            &mut post,
        );

        let probs: Vec<String> = post
            .into_iter()
            .rev()
            .map(|x| format!("{:.4}", x.2))
            .collect();

        assert_eq!(
            probs,
            ["0.5616", "0.4384", "0.8942", "0.1058", "0.9050", "0.0950"]
        );
    }
}
