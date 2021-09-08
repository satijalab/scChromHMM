# scChromHMM [![Rust](https://github.com/satijalab/scChromHMM/actions/workflows/rust.yml/badge.svg)](https://github.com/satijalab/scChromHMM/actions/workflows/rust.yml)

`scChromHMM` provides a suite of tools for the rapid processing of single-cell histone modification data to perform chromatin states analysis of the genome within each single-cell. It is an extention of bulk ChromHMM framework, which consumes the HMM model learned from ChromHMM and perform chromatin state analysis by running forward-backward algorithm for each single-cell.

# Input Data
scChromHMM primarily requires a group of four kind of files, which are defined as follows: 
* _fragment files_: Fragment files contains the information about the mapping location of a sequencing read fragments on the genome. The basic format is similar to as described by 10x, and it's primarily a BED file with an additional information of cellular barcode for each mapped fragment. example: `h3k27ac_fragments.tsv.gz`
    * **NOTE** the tabix index of the fragment files is also needed and can be generated using the command `tabix -f -p bed <fragment_file.gz>` for a block zipped (bgzip) fragment file. `example/h3k27ac_fragments.tsv.gz.tbi`
* _hmm_model_: A tsv file containing the information about the hmm model parameters. The default schema of this file is similar to the one generated by ChromHMM. `example/model_2.txt`.
* _anchors_: A tsv file with the list of anchors from the query data onto the reference data, along with their anchroring scores. `example/k27ac.txt`.
* _reference_cells_: A list of all the cellular barcodes (one per line) present in the reference dataset. `example/cells.txt`

# Compilation of the program
scChromHMM has been tested with stable release 1.52.1 of Rust, and the program can be compiled by using the command:

```{bash}
$ cargo build --release
```

# Running scChromHMM
Once compiled the scChromHMM program can be run to generate the posterior probability distribution across the hidden states using the command:
```{bash}
$ target/release/schrom hmm -f <fragment_files> -m <hmm_model> -a <anchor_files> -c <reference_cells> -t <number_of_threads> -o <output_folder>
```
**Note**: The order of fragment files should be the same as the anchor files. A toy example can be run using the data present in the example folder using the following command: (An extra flag `--onlyone` has been added to run the toy example on a subsequence of chromosome 1).
```
RUST_BACKTRACE=full RUST_LOG="trace" /usr/bin/time target/release/schrom hmm -f example/h3k27ac_fragments.tsv.gz example/h3k27me3_fragments.tsv.gz example/h3k4me1_fragments.tsv.gz -m example/model_2.txt -a example/k27ac.txt example/k27me3.txt example/k4me1.txt -c example/cells.txt -t 10 -o output --onlyone
```

# Transforming cell-wise probabilities into state-wise "short" representation

# Importing the output posterior probabilities into R

