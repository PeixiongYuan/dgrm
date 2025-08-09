
use crate::error::{DgrmError, Result};
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use fast_float::parse as fast_parse;
use std::collections::HashSet;
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
            let mut row_data = Vec::with_capacity(record.len().saturating_sub(1));
            for i in 1..record.len() {
                let value_str = record.get(i).unwrap_or("");
                // Fast path parse; fall back to standard parser for special values like NaN/inf
                let value: f64 = fast_parse(value_str)
                    .or_else(|_| value_str.parse::<f64>())
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
        let dosage_matrix = Array2::from_shape_vec((n_samples, n_variants), dosage_data)
            .map_err(|e| DgrmError::InvalidFormat(format!("Matrix creation failed: {}", e)))?;
        
        // Report missing values but do not impute here; let computation phase handle it
        let missing_count = dosage_matrix.iter().filter(|&&x| !x.is_finite() || x.is_nan()).count();
        if missing_count > 0 {
            log::warn!("Found {} missing or non-finite values; they will be imputed during computation", missing_count);
        }

        // Data quality information
        let mut min_val = f64::INFINITY;
        let mut max_val = f64::NEG_INFINITY;
        let mut sum_val = 0.0f64;
        let mut cnt_val: usize = 0;
        for &x in dosage_matrix.iter() {
            if x.is_finite() {
                if x < min_val { min_val = x; }
                if x > max_val { max_val = x; }
                sum_val += x;
                cnt_val += 1;
            }
        }
        let mean_val = if cnt_val > 0 { sum_val / cnt_val as f64 } else { f64::NAN };
        
        log::info!("Data range: [{:.4}, {:.4}]", min_val, max_val);
        log::info!("Mean dosage: {:.4}", mean_val);

        Ok(VntrData {
            sample_ids,
            dosage_matrix,
            n_variants,
        })
    }

    /// Filter samples by a set of sample IDs and return a new VntrData
    pub fn filter_by_sample_set(&self, keep: &HashSet<String>) -> Result<Self> {
        if keep.is_empty() {
            return Ok(self.clone());
        }

        // Collect indices to keep
        let mut indices: Vec<usize> = Vec::new();
        let mut new_sample_ids: Vec<String> = Vec::new();
        for (idx, sid) in self.sample_ids.iter().enumerate() {
            if keep.contains(sid) {
                indices.push(idx);
                new_sample_ids.push(sid.clone());
            }
        }

        if indices.is_empty() {
            return Err(DgrmError::InvalidFormat(
                "No samples from --keep file were found in the input data".to_string(),
            ));
        }

        // Build new dosage matrix with selected rows
        let n_variants = self.n_variants;
        let mut new_matrix = Array2::<f64>::zeros((indices.len(), n_variants));
        for (new_i, &old_i) in indices.iter().enumerate() {
            let row = self.dosage_matrix.row(old_i);
            new_matrix.row_mut(new_i).assign(&row);
        }

        log::info!(
            "Applied --keep filter: kept {} of {} samples",
            new_sample_ids.len(),
            self.sample_ids.len()
        );

        Ok(VntrData {
            sample_ids: new_sample_ids,
            dosage_matrix: new_matrix,
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

/// Read a keep file (two-column FID IID or single-column ID) and return a set of IDs (use IID if two columns)
pub fn read_keep_file<P: AsRef<Path>>(path: P) -> Result<HashSet<String>> {
    let path = path.as_ref();
    if !path.exists() {
        return Err(DgrmError::FileNotFound { path: path.display().to_string() });
    }

    let content = std::fs::read_to_string(path)?;
    let mut set = HashSet::new();
    for (ln, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() { continue; }
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() == 1 {
            set.insert(cols[0].to_string());
        } else if cols.len() >= 2 {
            // Use second column (IID) by convention
            set.insert(cols[1].to_string());
        } else {
            return Err(DgrmError::InvalidFormat(format!("Invalid line {} in keep file", ln + 1)));
        }
    }

    log::info!("Loaded {} sample IDs from --keep file: {}", set.len(), path.display());
    Ok(set)
}

