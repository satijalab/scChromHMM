#[macro_use]
extern crate log;

use clap::{App, Arg, SubCommand};
use std::error::Error;

mod config;
mod fragment;
mod hmm;
mod model;
mod quantify;
mod record;
mod transform;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("schrom")
        .version("0.1.0")
        .author("Avi Srivastava, Bingjie Zhang, Rahul Satija")
        .about("Generate summary stats for multimodal data.")
        .subcommand(
            SubCommand::with_name("hmm")
                .about("A subcommand to run hmm.")
                .arg(
                    Arg::with_name("fragments")
                        .long("fragments")
                        .short("f")
                        .takes_value(true)
                        .required(true)
                        .multiple(true)
                        .help("path to the fragment files."),
                )
                .arg(
                    Arg::with_name("anchors")
                        .long("anchors")
                        .short("a")
                        .takes_value(true)
                        .required(true)
                        .multiple(true)
                        .help("path to the anchors files. [Same order as fragments]"),
                )
                .arg(
                    Arg::with_name("threads")
                        .long("threads")
                        .short("t")
                        .takes_value(true)
                        .required(true)
                        .help("number of threads to use"),
                )
                .arg(
                    Arg::with_name("onlyone")
                        .help("quantify only chromosome one")
                        .long("onlyone")
                )
                .arg(
                    Arg::with_name("model")
                        .long("model")
                        .short("m")
                        .takes_value(true)
                        .required(true)
                        .help("path to the chromeHMM model.txt file"),
                )
                .arg(
                    Arg::with_name("output")
                        .long("output")
                        .short("o")
                        .takes_value(true)
                        .required(true)
                        .help("path to the output directory"),
                )
                .arg(
                    Arg::with_name("common_cells")
                        .long("common_cells")
                        .short("c")
                        .takes_value(true)
                        .required(true)
                        .help("path to the file with cellular barcodes of common assay"),
                ),
        )
        .subcommand(
            SubCommand::with_name("transform")
                .about("A subcommand to transform long form matrices to short.")
                .arg(
                    Arg::with_name("in_directory")
                        .long("in_directory")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the scChromHMM output directory"),
                )
                .arg(
                    Arg::with_name("out_directory")
                        .long("out_directory")
                        .short("o")
                        .takes_value(true)
                        .required(true)
                        .help("path to the scChromHMM output directory"),
                )
                .arg(
                    Arg::with_name("common_cells")
                        .long("common_cells")
                        .short("c")
                        .takes_value(true)
                        .required(true)
                        .help("path to the file with cellular barcodes of common assay."),
                ),
        )
        .get_matches();
    pretty_env_logger::init_timed();

    if let Some(sub_m) = matches.subcommand_matches("hmm") {
        hmm::callback(&sub_m)?
    }

    if let Some(sub_m) = matches.subcommand_matches("transform") {
        transform::callback(&sub_m)?
    }

    Ok(())
}
