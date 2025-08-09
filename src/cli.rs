
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "dgrm",
    version = "0.1.0",
    about = "Parallel VNTR GRM calculator implemented in Rust",
    long_about = r#"
VNTR GRM Calculator - Rust Implementation

This tool calculates Genetic Relationship Matrix (GRM) from VNTR (Variable Number Tandem Repeat) 
scaled dosage data and outputs in GCTA-compatible format.

Input file format:
  - No header
  - First column: sample IDs
  - Other columns: scaled dosage values for each variant
  - Rows: samples, Columns: variants

Output: GCTA compatible GRM files
  - <prefix>.grm.N.bin: Number of variants (binary)
  - <prefix>.grm.bin: GRM values (binary)
  - <prefix>.grm.id: Sample IDs (text)
  - <prefix>.grm.summary: Summary statistics (text)

Example:
  dgrm --input chr1_vntr_scaled_dos.txt --output chr1_grm --threads 8
"#
)]
pub struct Args {
    /// Input VNTR scaled dosage file
    #[arg(
        short = 'i',
        long = "input",
        value_name = "FILE",
        help = "Path to input VNTR scaled dosage file"
    )]
    pub input: PathBuf,

    /// Output prefix for GCTA format files
    #[arg(
        short = 'o',
        long = "output",
        value_name = "PREFIX",
        help = "Output prefix for GCTA format files"
    )]
    pub output: String,

    /// Number of threads for parallel computation
    #[arg(
        short = 't',
        long = "threads",
        value_name = "N",
        default_value = "1",
        help = "Number of threads for parallel computation"
    )]
    pub threads: usize,

    /// Keep only samples listed in this file (two columns: FID IID or single column)
    #[arg(
        long = "keep",
        value_name = "KEEP_FILE",
        help = "Path to a file containing sample IDs to keep (format: FID IID or single ID per line)"
    )]
    pub keep: Option<PathBuf>,

    

}