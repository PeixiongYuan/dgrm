
use crate::error::{DgrmError, Result};
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use std::fs::File;
use std::path::Path;
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct VntrData {
    /// Sample IDs
    sample_ids: Vec<String>,
    /// Dosage matrix: rows = samples, columns = variants
    dosage_matrix: Array2<f64>,
    /// Number of variants (more efficient than storing IDs)
    n_variants: usize,
}

impl VntrData {
    /// Load VNTR data from scaled dosage file
    pub fn load<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        let file_path = file_path.as_ref();
        
        if !file_path.exists() {
            return Err(DgrmError::FileNotFound {
                path: file_path.display().to_string(),
            });
        }

        // Check file size
        let file_size = std::fs::metadata(file_path)?.len();
        let file_size_mb = file_size as f64 / (1024.0 * 1024.0);
        log::info!("File size: {:.2} MB", file_size_mb);

        let start_time = Instant::now();
        
        // Open file
        let file = File::open(file_path)?;
        let mut reader = ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(file);

        // Read all records
        let mut sample_ids = Vec::new();
        let mut dosage_data = Vec::new();
        let mut n_variants = 0;

        // Setup progress bar
        let pb = ProgressBar::new_spinner();
        pb.set_style(ProgressStyle::default_spinner()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ")
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap());
        pb.set_message("Reading data...");

        for (line_num, result) in reader.records().enumerate() {
            let record = result?;
            
            if record.is_empty() {
                continue;
            }

            // First column is sample ID
            let sample_id = record.get(0)
                .ok_or_else(|| DgrmError::InvalidFormat(
                    format!("Missing sample ID on line {}", line_num + 1)
                ))?;
            sample_ids.push(sample_id.to_string());

            // Remaining columns are dosage values
            let mut row_data = Vec::new();
            for i in 1..record.len() {
                let value_str = record.get(i).unwrap_or("0");
                let value: f64 = value_str.parse()
                    .map_err(|_| DgrmError::InvalidFormat(
                        format!("Invalid numeric value '{}' on line {}, column {}", 
                               value_str, line_num + 1, i + 1)
                    ))?;
                row_data.push(value);
            }

            // Check consistency of variant count
            if n_variants == 0 {
                n_variants = row_data.len();
            } else if row_data.len() != n_variants {
                return Err(DgrmError::DimensionMismatch {
                    expected: n_variants,
                    actual: row_data.len(),
                });
            }

            dosage_data.extend(row_data);

            if line_num % 1000 == 0 {
                pb.set_message(format!("Read {} samples", line_num + 1));
                pb.tick();
            }
        }

        pb.finish_with_message("Data reading completed");

        let n_samples = sample_ids.len();
        let read_time = start_time.elapsed();
        
        log::info!("Data reading completed in {:.2} seconds", read_time.as_secs_f64());
        log::info!("Data dimensions: {} samples x {} variants", n_samples, n_variants);

        if n_samples == 0 || n_variants == 0 {
            return Err(DgrmError::InvalidFormat(
                "No valid data found in input file".to_string()
            ));
        }

        // Create dosage matrix
        let mut dosage_matrix = Array2::from_shape_vec((n_samples, n_variants), dosage_data)
            .map_err(|e| DgrmError::InvalidFormat(format!("Matrix creation failed: {}", e)))?;

        // Handle missing values (NaN) - replace with 0
        let missing_count = dosage_matrix.iter().filter(|&&x| x.is_nan()).count();
        if missing_count > 0 {
            log::warn!("Found {} missing values, replacing with 0", missing_count);
            dosage_matrix.mapv_inplace(|x| if x.is_nan() { 0.0 } else { x });
        }

        // Data quality information
        let min_val = dosage_matrix.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_val = dosage_matrix.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let mean_val = dosage_matrix.mean().unwrap_or(0.0);
        
        log::info!("Data range: [{:.4}, {:.4}]", min_val, max_val);
        log::info!("Mean dosage: {:.4}", mean_val);

        Ok(VntrData {
            sample_ids,
            dosage_matrix,
            n_variants,
        })
    }

    /// Get number of samples
    pub fn n_samples(&self) -> usize {
        self.sample_ids.len()
    }

    /// Get number of variants
    pub fn n_variants(&self) -> usize {
        self.n_variants
    }

    /// Get sample IDs
    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    /// Get dosage matrix
    pub fn dosage_matrix(&self) -> &Array2<f64> {
        &self.dosage_matrix
    }
}

