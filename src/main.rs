use clap::{Command, Arg, ArgMatches};

mod parallel_toy;
mod parallel_toy_ipc;

fn main() {
    let matches = Command::new("Toy Parallel BAM Converter")
        .version("0.1.0")
        .about("Test parallel processing for BAM to Parquet/Arrow IPC conversion")
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .value_name("FILE")
            .help("Output file path (.parquet or .arrow)")
            .required(true))
        .arg(Arg::new("format")
            .short('f')
            .long("format")
            .value_name("FORMAT")
            .help("Output format: parquet or ipc")
            .default_value("parquet"))
        .arg(Arg::new("batches")
            .short('b')
            .long("batches")
            .value_name("NUMBER")
            .help("Number of batches to process")
            .default_value("10"))
        .arg(Arg::new("records")
            .short('r')
            .long("records")
            .value_name("NUMBER")
            .help("Records per batch")
            .default_value("10000"))
        .arg(Arg::new("threads")
            .short('t')
            .long("threads")
            .value_name("NUMBER")
            .help("Number of threads")
            .default_value("4"))
        .arg(Arg::new("no-sequence")
            .long("no-sequence")
            .help("Exclude sequence data")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("no-quality")
            .long("no-quality")
            .help("Exclude quality scores")
            .action(clap::ArgAction::SetTrue))
        .get_matches();

    let output_path = matches.get_one::<String>("output").unwrap();
    let format = matches.get_one::<String>("format").unwrap();
    let num_batches: usize = matches.get_one::<String>("batches").unwrap().parse().expect("Invalid number of batches");
    let records_per_batch: usize = matches.get_one::<String>("records").unwrap().parse().expect("Invalid records per batch");
    let num_threads: usize = matches.get_one::<String>("threads").unwrap().parse().expect("Invalid number of threads");
    let include_sequence = !matches.get_flag("no-sequence");
    let include_quality = !matches.get_flag("no-quality");

    println!("Running parallel toy conversion with {} format...", format);
    
    let result = match format.as_str() {
        "ipc" | "arrow" => {
            parallel_toy_ipc::parallel_toy_ipc_conversion(
                output_path,
                num_batches,
                records_per_batch,
                num_threads,
                include_sequence,
                include_quality,
            )
        }
        "parquet" => {
            parallel_toy::parallel_toy_conversion(
                output_path,
                num_batches,
                records_per_batch,
                num_threads,
                include_sequence,
                include_quality,
            )
        }
        _ => {
            eprintln!("Unsupported format: {}. Use 'parquet' or 'ipc'", format);
            std::process::exit(1);
        }
    };
    
    match result {
        Ok(()) => println!("Success!"),
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}
