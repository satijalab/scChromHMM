# scChromHMM

`scChromHMM` provides a suite of tools for the rapid processing of single-cell histone modification data to perform chromatin states analysis of the genome within each single-cell. It is an extention of bulk ChromHMM framework, which consumes the HMM model learned from ChromHMM and perform chromatin state analysis by running forward-backward algorithm for each single-cell.

# Building from source
scChromHMM has been build with stable release 1.52.1 of Rust, and the program can be compiled by using the command:

```{bash}
$ cargo build --release
```

