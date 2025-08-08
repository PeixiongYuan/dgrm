
use anyhow::Result;
use clap::Parser;

mod cli;
mod data_io;
mod grm_calculator;
mod gcta_writer;
mod error;

use grm_calculator::MafFilter;

use cli::Args;
use data_io::VntrData;
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
    let vntr_data = VntrData::load(&args.input)?;
    log::info!("Loaded {} samples and {} variants", 
               vntr_data.n_samples(), vntr_data.n_variants());
    
    // Prepare MAF filter
    let maf_filter = if args.maf.is_some() || args.max_maf.is_some() {
        Some(MafFilter {
            min_maf: args.maf,
            max_maf: args.max_maf,
        })
    } else {
        None
    };

    if let Some(ref filter) = maf_filter {
        if let Some(min_maf) = filter.min_maf {
            log::info!("Applying minimum MAF filter: {:.4}", min_maf);
        }
        if let Some(max_maf) = filter.max_maf {
            log::info!("Applying maximum MAF filter: {:.4}", max_maf);
        }
    } else {
        log::info!("No MAF filtering applied");
    }

    // Calculate GRM
    log::info!("Calculating GRM matrix...");
    let calculator = GrmCalculator::new(args.threads);
    let grm_matrix = calculator.calculate_grm(&vntr_data, maf_filter)?;
    
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