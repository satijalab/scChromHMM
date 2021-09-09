# scChromHMM [![Rust](https://github.com/satijalab/scChromHMM/actions/workflows/rust.yml/badge.svg)](https://github.com/satijalab/scChromHMM/actions/workflows/rust.yml)

`scChromHMM` provides a suite of tools for the rapid processing of single-cell histone modification data to perform chromatin states analysis of the genome within each single-cell. It is an extention of bulk ChromHMM framework, which consumes the HMM model learned from ChromHMM and perform chromatin state analysis by running forward-backward algorithm for each single-cell.

# Input Data
scChromHMM primarily requires a group of four kind of files, which are defined as follows: 
* _fragment files_: Fragment files contains the information about the mapping location of a sequencing read fragments on the genome. The basic format is similar to as described by 10x, and it's primarily a BED file with an additional information of cellular barcode for each mapped fragment. toy example: `h3k27ac_fragments.tsv.gz`
    * **NOTE** the tabix index of the fragment files is also needed and can be generated using the command `tabix -f -p bed <fragment_file.gz>` for a block zipped (bgzip) fragment file. toy example: `example/h3k27ac_fragments.tsv.gz.tbi`
* _hmm_model_: A tsv file containing the information about the hmm model parameters. The default schema of this file is similar to the one generated by ChromHMM. toy example: `example/model_2.txt`.
* _anchors_: A tsv file with the list of anchors from the query data onto the reference data, along with their anchroring scores. toy example:`example/k27ac.txt`.
* _reference_cells_: A list of all the cellular barcodes (one per line) present in the reference dataset. toy example:`example/cells.txt`

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

# State-wise "short" representation
The `hmm` subcommand of the scChromHMM tool generates cell-wise posterior probabilities for every reference cell across the genome. The probabilities are stored for each cell in a binary format i.e. 200bp region by state matrix with integer values in range [0-100]. toy example: `output/chr1/L1_CCTCTAGTCGCTAAAC.bin`. Based on the number of reference cells, size of the output posterior probabilites can grow significantly; and some downstream analyses are faster to work with region by cells matrix (for each state) instead of region by state (for each cell) matrices. Hence, scChromHMM subcommand `transform` can be used to convert the data into the "short" representation of region by cell. The command to do that is as follows:
```{bash}
$ target/release/schrom transform -c <reference_cells> -i <input_folder> -o <output_folder>
```
The toy example can be run using the following command. **NOTE** An extra flag `--onlyone` has been added to run the toy example on a subsequence of chromosome 1.
```bash
$ mkdir short_output
$ RUST_BACKTRACE=full RUST_LOG="trace" /usr/bin/time target/release/schrom transform -c example/cells.txt -i output -o short_output --onlyone
```

# Importing the posterior probabilities into R
The chromatin state wise, region by cells posterior probabilities of the toy example can be imported into an R environment using YYYY package as follows:
```{R}
library(YYY)
mat <- get_indus_state("short_output/chr1/1.bin", "chr1", "short_output/chr1/cells.txt")
dim(mat)
# [1] 5001 7201
```

