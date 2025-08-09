
use anyhow::Result;
use clap::Parser;

mod cli;
mod data_io;
mod grm_calculator;
mod gcta_writer;
mod error;


use cli::Args;
use data_io::{VntrData, read_keep_file};
use grm_calculator::GrmCalculator;
use gcta_writer::GctaWriter;

fn main() -> Result<()> {
    env_logger::init();
    
    let args = Args::parse();
    
    log::info!("=== VNTR GRM Calculator (Rust Implementation) ===");
    log::info!("Input file: {}", args.input.display());
    log::info!("Output prefix: {}", args.output);
    log::info!("Threads: {}", args.threads);
    
    // Set up rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .expect("Failed to initialize thread pool");
    
    // Load VNTR data
    log::info!("Loading VNTR data...");
    let mut vntr_data = VntrData::load(&args.input)?;
    log::info!("Loaded {} samples and {} variants", 
               vntr_data.n_samples(), vntr_data.n_variants());

    // Apply --keep filter if provided
    if let Some(keep_path) = &args.keep {
        let keep_set = read_keep_file(keep_path)?;
        vntr_data = vntr_data.filter_by_sample_set(&keep_set)?;
        log::info!(
            "After --keep filtering: {} samples remain",
            vntr_data.n_samples()
        );
    }
    
    // No MAF filtering
    log::info!("No MAF filtering (parameter removed)");

    // Calculate GRM
    log::info!("Calculating GRM matrix...");
    let calculator = GrmCalculator::new(args.threads);
    let grm_matrix = calculator.calculate_grm(&vntr_data)?;
    
    // Write GCTA format files
    log::info!("Writing GCTA format files...");
    let writer = GctaWriter::new(&args.output);
    writer.write_grm(&grm_matrix, &vntr_data)?;
    
    log::info!("=== SUCCESS ===");
    log::info!("VNTR GRM calculation completed successfully!");
    log::info!("Total samples processed: {}", vntr_data.n_samples());
    log::info!("Total variants processed: {}", vntr_data.n_variants());
    
    Ok(())
}