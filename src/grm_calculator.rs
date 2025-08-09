
use crate::data_io::VntrData;
use crate::error::Result;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use ndarray::s;
use rayon::prelude::*;
use std::sync::Arc;
use std::time::Instant;

/// MAF filtering options
#[derive(Debug, Clone)]
pub struct MafFilter {
    pub min_maf: Option<f64>,
    pub max_maf: Option<f64>,
}

impl Default for MafFilter {
    fn default() -> Self {
        Self {
            min_maf: None,
            max_maf: None,
        }
    }
}

pub struct GrmCalculator {
    threads: usize,
}

#[derive(Debug, Clone)]
pub struct GrmMatrix {
    /// Correlation matrix (symmetric)
    matrix: Array2<f64>,
    /// Sample IDs
    sample_ids: Vec<String>,
    /// Number of variants used in calculation
    n_variants: usize,
}

impl GrmCalculator {
    pub fn new(threads: usize) -> Self {
        Self { threads }
    }

    /// Calculate GRM matrix using standard genetic relationship matrix formula
    /// GRM_ij = (1/M) × Σ_k [(x_ik - 2p_k)(x_jk - 2p_k)] / [2p_k(1-p_k)]
    pub fn calculate_grm(&self, vntr_data: &VntrData, maf_filter: Option<MafFilter>) -> Result<GrmMatrix> {
        self.calculate_grm_with_filter(vntr_data, maf_filter.unwrap_or_default())
    }

    /// Calculate GRM with MAF filtering
    fn calculate_grm_with_filter(&self, vntr_data: &VntrData, maf_filter: MafFilter) -> Result<GrmMatrix> {
        let start_time = Instant::now();
        let n_samples = vntr_data.n_samples();
        let n_variants = vntr_data.n_variants();
        
        log::info!("Starting GRM calculation for {} samples with {} variants", 
                   n_samples, n_variants);
        log::info!("Using {} threads", self.threads);

        // Comprehensive memory estimation
        let grm_elements = (n_samples * (n_samples + 1)) / 2;
        let full_matrix_elements = n_samples * n_samples;
        
        // Memory breakdown:
        let grm_matrix_mb = (full_matrix_elements * 8) / (1024 * 1024); // Full symmetric matrix
        let dosage_matrix_mb = (n_samples * n_variants * 8) / (1024 * 1024); // Input data
        let sample_stats_mb = (n_samples * 2 * 8) / (1024 * 1024); // Stats cache
        let correlation_temp_mb = (grm_elements * 16) / (1024 * 1024); // Temp correlation vector
        let overhead_mb = 500; // System overhead, thread stacks, etc.
        
        let total_estimated_mb = grm_matrix_mb + dosage_matrix_mb + sample_stats_mb + correlation_temp_mb + overhead_mb;
        
        log::info!("Detailed memory estimation:");
        log::info!("  GRM matrix ({} x {}): {} MB", n_samples, n_samples, grm_matrix_mb);
        log::info!("  Dosage matrix ({} x {}): {} MB", n_samples, n_variants, dosage_matrix_mb);
        log::info!("  Temporary correlation data: {} MB", correlation_temp_mb);
        log::info!("  Sample statistics cache: {} MB", sample_stats_mb);
        log::info!("  System overhead: {} MB", overhead_mb);
        log::info!("  Total estimated memory: {} MB ({:.1} GB)", total_estimated_mb, total_estimated_mb as f64 / 1024.0);
        
        // Ultra-aggressive memory safety for very large datasets
        // Experience shows that actual memory usage can be 2-3x higher than estimates
        if n_samples > 30000 || total_estimated_mb > 20 * 1024 { // > 30K samples or > 20GB
            log::warn!("⚠️  Large dataset detected: {} samples, estimated {:.1} GB memory", n_samples, total_estimated_mb as f64 / 1024.0);
            log::warn!("⚠️  Actual memory usage often exceeds estimates by 2-3x in parallel computation");
            log::info!("🔄 Automatically using memory-efficient block computation for safety");
            
            // Calculate conservative block size to keep peak memory under 40GB
            let target_memory_gb = 40.0;
            let target_memory_mb = (target_memory_gb * 1024.0) as usize;
            let available_for_grm = target_memory_mb - dosage_matrix_mb - overhead_mb;
            let max_grm_elements = available_for_grm * 1024 * 1024 / 8; // Convert to elements
            let max_samples_per_block = ((max_grm_elements as f64).sqrt() as usize).max(1024).min(4096); // 1K-4K range
            
            log::info!("Block size: {} samples per block (target memory: {:.1} GB)", max_samples_per_block, target_memory_gb);
            return self.calculate_grm_blockwise(vntr_data, max_samples_per_block, maf_filter);
        } else if total_estimated_mb > 15 * 1024 { // > 15GB, warn but continue  
            log::warn!("⚠️  MODERATE MEMORY WARNING: Estimated usage {:.1} GB", total_estimated_mb as f64 / 1024.0);
            log::warn!("⚠️  Actual usage may be higher - consider monitoring system memory");
        }

        // Create progress bar
        let pb = ProgressBar::new(grm_elements as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
            .unwrap()
            .progress_chars("#>-"));
        pb.set_message("Calculating correlations");

        let dosage_matrix = vntr_data.dosage_matrix();
        
        // Step 1: Calculate allele frequencies for each variant
        log::info!("Calculating allele frequencies for {} variants...", n_variants);
        let allele_freqs: Vec<f64> = (0..n_variants)
            .into_par_iter()
            .map(|variant_idx| {
                // Calculate allele frequency for this variant across all samples
                let mut sum = 0.0;
                let mut valid_count = 0;
                
                for sample_idx in 0..n_samples {
                    let dosage = dosage_matrix[[sample_idx, variant_idx]];
                    if !dosage.is_nan() && dosage.is_finite() {
                        sum += dosage;
                        valid_count += 1;
                    }
                }
                
                if valid_count > 0 {
                    // Allele frequency = mean dosage / 2 (since dosage is 0, 1, or 2)
                    (sum / valid_count as f64) / 2.0
                } else {
                    0.5 // Default frequency if no valid data
                }
            })
            .collect();
        
        // Step 2: Filter variants based on MAF thresholds
        let min_maf = maf_filter.min_maf.unwrap_or(0.0);
        let max_maf = maf_filter.max_maf.unwrap_or(1.0);
        
        let valid_variants: Vec<usize> = allele_freqs
            .iter()
            .enumerate()
            .filter(|(_, &freq)| {
                // Calculate MAF (minor allele frequency)
                let maf = freq.min(1.0 - freq);
                maf >= min_maf && maf <= max_maf
            })
            .map(|(idx, _)| idx)
            .collect();
        
        let n_valid_variants = valid_variants.len();
        if maf_filter.min_maf.is_some() || maf_filter.max_maf.is_some() {
            log::info!("Using {} valid variants (MAF: {:.4} - {:.4}) out of {} total variants", 
                       n_valid_variants, min_maf, max_maf, n_variants);
        } else {
            log::info!("Using {} variants (no MAF filtering) out of {} total variants", 
                       n_valid_variants, n_variants);
        }
        
        if n_valid_variants == 0 {
            return Err(crate::error::DgrmError::InvalidFormat(
                "No valid variants found (all variants are monomorphic)".to_string()
            ));
        }

        // Build imputed data matrix for valid variants only
        log::info!("Computing GRM via standardized matrix multiplication...");
        pb.set_message("Building imputed matrix Z");
        let mut z = Array2::<f64>::zeros((n_samples, n_valid_variants));
        for (new_idx, &variant_idx) in valid_variants.iter().enumerate() {
            let mean_val = allele_freqs[variant_idx] * 2.0;
            for sample_idx in 0..n_samples {
                let value = dosage_matrix[[sample_idx, variant_idx]];
                z[[sample_idx, new_idx]] = if value.is_finite() && !value.is_nan() {
                    value
                } else {
                    mean_val
                };
            }
        }

        // Standardize rows (samples): subtract mean, divide by std
        pb.set_message("Standardizing rows");
        let stats: Vec<(f64, f64)> = (0..n_samples)
            .into_par_iter()
            .map(|i| {
                let row = z.row(i);
                let mean = row.sum() / n_valid_variants as f64;
                let var = row
                    .iter()
                    .map(|&x| {
                        let d = x - mean;
                        d * d
                    })
                    .sum::<f64>()
                    / (n_valid_variants - 1) as f64;
                (mean, var.sqrt())
            })
            .collect();

        for i in 0..n_samples {
            let (mean, std_raw) = stats[i];
            let std = std_raw.max(f64::EPSILON);
            let mut row = z.row_mut(i);
            for x in row.iter_mut() {
                *x = (*x - mean) / std;
            }
        }

        // Compute GRM = Z * Z^T / (m-1) using upper-triangular blocked parallel matmul
        pb.set_message("Matrix multiply Z * Z^T (blocked upper-triangular)");
        let z_t = z.t().to_owned();
        let chunk: usize = std::env::var("DGRM_MM_CHUNK").ok()
            .and_then(|v| v.parse::<usize>().ok()).unwrap_or(512);

        // Prepare block coordinates for upper triangle (including diagonal)
        let mut coords: Vec<(usize, usize, usize, usize)> = Vec::new();
        let mut rs = 0usize;
        while rs < n_samples {
            let re = (rs + chunk).min(n_samples);
            let mut cs = rs;
            while cs < n_samples {
                let ce = (cs + chunk).min(n_samples);
                coords.push((rs, re, cs, ce));
                cs += chunk;
            }
            rs += chunk;
        }

        let parts: Vec<((usize, usize, usize, usize), ndarray::Array2<f64>)> = coords
            .into_par_iter()
            .map(|(rs, re, cs, ce)| {
                let a = z.slice(s![rs..re, ..]);
                let b = z_t.slice(s![.., cs..ce]);
                let sub = a.dot(&b);
                ((rs, re, cs, ce), sub)
            })
            .collect();

        let mut grm_matrix = Array2::<f64>::zeros((n_samples, n_samples));
        for ((rs, re, cs, ce), sub) in parts.into_iter() {
            let rows = re - rs;
            let cols = ce - cs;
            // write upper block
            for i in 0..rows {
                for j in 0..cols {
                    grm_matrix[[rs + i, cs + j]] = sub[[i, j]];
                }
            }
            // mirror to lower block if off-diagonal
            if cs != rs {
                for i in 0..rows {
                    for j in 0..cols {
                        grm_matrix[[cs + j, rs + i]] = sub[[i, j]];
                    }
                }
            }
        }
        let scale = 1.0 / (n_valid_variants.saturating_sub(1)) as f64;
        grm_matrix.mapv_inplace(|x| x * scale);
        for i in 0..n_samples { grm_matrix[[i, i]] = 1.0; }

        pb.finish_with_message("GRM calculation completed");

        let calc_time = start_time.elapsed();
        log::info!("GRM calculation completed in {:.2} seconds", calc_time.as_secs_f64());
        
        // Calculate basic statistics more efficiently
        let upper_triangular = self.extract_upper_triangular(&grm_matrix);
        let n_elements = upper_triangular.len();
        let sum = upper_triangular.iter().sum::<f64>();
        let mean_kinship = sum / n_elements as f64;
        
        let variance = upper_triangular.iter()
            .map(|&x| (x - mean_kinship).powi(2))
            .sum::<f64>() / n_elements as f64;
        let sd_kinship = variance.sqrt();
        
        let min_kinship = upper_triangular.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_kinship = upper_triangular.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        log::info!("GRM Statistics:");
        log::info!("  Mean kinship: {:.6}", mean_kinship);
        log::info!("  SD kinship: {:.6}", sd_kinship);
        log::info!("  Min kinship: {:.6}", min_kinship);
        log::info!("  Max kinship: {:.6}", max_kinship);

        Ok(GrmMatrix {
            matrix: grm_matrix,
            sample_ids: vntr_data.sample_ids().to_vec(),
            n_variants: n_valid_variants, // Use valid variant count
        })
    }




    /// Extract upper triangular values for statistics calculation
    /// Uses column-major order to match GCTA/R expectations
    fn extract_upper_triangular(&self, matrix: &Array2<f64>) -> Vec<f64> {
        let n = matrix.nrows();
        let mut values = Vec::with_capacity((n * (n + 1)) / 2);
        
        // Column-major order: for each column j, add elements from row j to n-1
        for j in 0..n {
            for i in j..n {
                values.push(matrix[[i, j]]);
            }
        }
        
        values
    }

    /// Memory-efficient block-wise GRM calculation for large datasets
    pub fn calculate_grm_blockwise(&self, vntr_data: &VntrData, block_size: usize, maf_filter: MafFilter) -> Result<GrmMatrix> {
        let start_time = Instant::now();
        let n_samples = vntr_data.n_samples();
        let n_variants = vntr_data.n_variants();
        
        log::info!("Starting memory-efficient block-wise GRM calculation");
        log::info!("Total samples: {}, Block size: {}", n_samples, block_size);
        
        // Calculate number of blocks
        let n_blocks = (n_samples + block_size - 1) / block_size;
        let total_block_pairs = (n_blocks * (n_blocks + 1)) / 2;
        
        log::info!("Using {} blocks, {} block pairs to compute", n_blocks, total_block_pairs);
        
        // Initialize result matrix (this is the main memory consumer)
        let mut grm_matrix = Array2::<f64>::zeros((n_samples, n_samples));
        let dosage_matrix = vntr_data.dosage_matrix();
        
        // Calculate allele frequencies and filter variants (same as main algorithm)
        log::info!("Calculating allele frequencies for {} variants...", n_variants);
        let allele_freqs: Vec<f64> = (0..n_variants)
            .into_par_iter()
            .map(|variant_idx| {
                let mut sum = 0.0;
                let mut valid_count = 0;
                
                for sample_idx in 0..n_samples {
                    let dosage = dosage_matrix[[sample_idx, variant_idx]];
                    if !dosage.is_nan() && dosage.is_finite() {
                        sum += dosage;
                        valid_count += 1;
                    }
                }
                
                if valid_count > 0 {
                    (sum / valid_count as f64) / 2.0
                } else {
                    0.5
                }
            })
            .collect();
        
        // Filter variants based on MAF thresholds
        let min_maf = maf_filter.min_maf.unwrap_or(0.0);
        let max_maf = maf_filter.max_maf.unwrap_or(1.0);
        
        let valid_variants: Vec<usize> = allele_freqs
            .iter()
            .enumerate()
            .filter(|(_, &freq)| {
                // Calculate MAF (minor allele frequency)
                let maf = freq.min(1.0 - freq);
                maf >= min_maf && maf <= max_maf
            })
            .map(|(idx, _)| idx)
            .collect();
        
        let n_valid_variants = valid_variants.len();
        if maf_filter.min_maf.is_some() || maf_filter.max_maf.is_some() {
            log::info!("Using {} valid variants (MAF: {:.4} - {:.4}) out of {} total variants", 
                       n_valid_variants, min_maf, max_maf, n_variants);
        } else {
            log::info!("Using {} variants (no MAF filtering) out of {} total variants", 
                       n_valid_variants, n_variants);
        }
        
        if n_valid_variants == 0 {
            return Err(crate::error::DgrmError::InvalidFormat(
                "No valid variants found (all variants are monomorphic)".to_string()
            ));
        }
        
        let valid_variants = Arc::new(valid_variants);
        let _allele_freqs = Arc::new(allele_freqs);
        
        // Progress tracking
        let pb = ProgressBar::new(total_block_pairs as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} blocks ({percent}%) {msg}")
            .unwrap()
            .progress_chars("#>-"));
        pb.set_message("Processing blocks");

        // Process blocks with aggressive memory management
        for block_i in 0..n_blocks {
            for block_j in block_i..n_blocks {
                let start_i = block_i * block_size;
                let end_i = (start_i + block_size).min(n_samples);
                let start_j = block_j * block_size;
                let end_j = (start_j + block_size).min(n_samples);

                // Process smaller sub-blocks within each block to reduce peak memory
                let sub_block_size = 512; // Process in 512-sample chunks
                for sub_i_start in (start_i..end_i).step_by(sub_block_size) {
                    let sub_i_end = (sub_i_start + sub_block_size).min(end_i);
                    
                    // Collect GRM values for this sub-block
                    let sub_block_grm_values: Vec<((usize, usize), f64)> = (sub_i_start..sub_i_end)
                        .into_par_iter()
                        .flat_map(|i| {
                            let start_j_inner = if block_i == block_j { i.max(start_j) } else { start_j };
                            let variants_ref = valid_variants.clone();
                            let afreqs = _allele_freqs.clone();
                            (start_j_inner..end_j).into_par_iter().map(move |j| {
                                // Calculate Pearson correlation for all pairs (including diagonal)
                                let mut sum_i = 0.0;
                                let mut sum_j = 0.0;
                                let mut sum_ij = 0.0;
                                let mut sum_i2 = 0.0;
                                let mut sum_j2 = 0.0;
                                let mut count = 0;
                                
                                for &variant_idx in variants_ref.iter() {
                                    let mean_val = afreqs[variant_idx] * 2.0;
                                    let mut x_i = dosage_matrix[[i, variant_idx]];
                                    let mut x_j = dosage_matrix[[j, variant_idx]];
                                    if !x_i.is_finite() || x_i.is_nan() { x_i = mean_val; }
                                    if !x_j.is_finite() || x_j.is_nan() { x_j = mean_val; }
                                    sum_i += x_i;
                                    sum_j += x_j;
                                    sum_ij += x_i * x_j;
                                    sum_i2 += x_i * x_i;
                                    sum_j2 += x_j * x_j;
                                    count += 1;
                                }
                                
                                let grm_value = if count > 1 {
                                    let n = count as f64;
                                    let numerator = n * sum_ij - sum_i * sum_j;
                                    let denom_i = n * sum_i2 - sum_i * sum_i;
                                    let denom_j = n * sum_j2 - sum_j * sum_j;
                                    
                                    if denom_i > 0.0 && denom_j > 0.0 {
                                        numerator / (denom_i * denom_j).sqrt()
                                    } else {
                                        // Handle zero variance case
                                        if i == j { 1.0 } else { 0.0 }
                                    }
                                } else {
                                    // Not enough data
                                    if i == j { 1.0 } else { 0.0 }
                                };
                                ((i, j), grm_value)
                            })
                        })
                        .collect();

                    // Immediately write results and drop temporary data
                    for ((i, j), grm_value) in sub_block_grm_values {
                        grm_matrix[[i, j]] = grm_value;
                        if i != j {
                            grm_matrix[[j, i]] = grm_value; // Symmetric
                        }
                    }
                    // sub_block_grm_values is dropped here, freeing memory
                }
                
                pb.inc(1);
                
                // Progress reporting and memory cleanup hints
                if (block_i * n_blocks + block_j) % 5 == 0 {
                    let progress_pct = (block_i * n_blocks + block_j + 1) as f64 / total_block_pairs as f64 * 100.0;
                    log::info!("Processed block pair ({}, {}) - {:.1}% complete", block_i, block_j, progress_pct);
                    
                    // Suggest garbage collection more frequently for large datasets
                    if n_samples > 40000 {
                        std::hint::black_box(&grm_matrix); // Prevent over-optimization
                    }
                }
            }
        }

        pb.finish_with_message("Block-wise calculation completed");
        
        let calc_time = start_time.elapsed();
        log::info!("Block-wise GRM calculation completed in {:.2} seconds", calc_time.as_secs_f64());

        // Calculate statistics efficiently
        let upper_triangular = self.extract_upper_triangular(&grm_matrix);
        let n_elements = upper_triangular.len();
        let sum = upper_triangular.iter().sum::<f64>();
        let mean_kinship = sum / n_elements as f64;
        
        let variance = upper_triangular.iter()
            .map(|&x| (x - mean_kinship).powi(2))
            .sum::<f64>() / n_elements as f64;
        let sd_kinship = variance.sqrt();
        
        let min_kinship = upper_triangular.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_kinship = upper_triangular.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        log::info!("GRM Statistics:");
        log::info!("  Mean kinship: {:.6}", mean_kinship);
        log::info!("  SD kinship: {:.6}", sd_kinship);
        log::info!("  Min kinship: {:.6}", min_kinship);
        log::info!("  Max kinship: {:.6}", max_kinship);

        Ok(GrmMatrix {
            matrix: grm_matrix,
            sample_ids: vntr_data.sample_ids().to_vec(),
            n_variants: n_valid_variants, // Use valid variant count
        })
    }


}

impl GrmMatrix {
    /// Get sample IDs
    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    /// Get number of variants
    #[allow(dead_code)]
    pub fn n_variants(&self) -> usize {
        self.n_variants
    }

    /// Get number of samples
    pub fn n_samples(&self) -> usize {
        self.sample_ids.len()
    }


    /// Get matrix value at (i, j)
    pub fn value(&self, i: usize, j: usize) -> f64 {
        self.matrix[[i, j]]
    }

    /// Get upper triangular values (including diagonal) for GCTA format
    /// Uses column-major order to match GCTA dspMatrix format (uplo='U')
    pub fn upper_triangular_values(&self) -> Vec<f64> {
        let n = self.matrix.nrows();
        let mut values = Vec::with_capacity((n * (n + 1)) / 2);
        
        // GCTA uses column-major order for upper triangular: 
        // for each column j, add elements from row 0 to j (upper triangle + diagonal)
        for j in 0..n {
            for i in 0..=j {
                values.push(self.matrix[[i, j]]);
            }
        }
        
        values
    }

    /// Calculate basic statistics
    pub fn calculate_stats(&self) -> GrmStats {
        let values = self.upper_triangular_values();
        let n_elements = values.len();
        
        let mean = values.iter().sum::<f64>() / n_elements as f64;
        let variance = values.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>() / n_elements as f64;
        let std_dev = variance.sqrt();
        
        let min_val = values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_val = values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        
        GrmStats {
            n_variants: self.n_variants,
            n_samples: self.n_samples(),
            n_grm_elements: n_elements,
            mean_kinship: mean,
            sd_kinship: std_dev,
            min_kinship: min_val,
            max_kinship: max_val,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GrmStats {
    pub n_variants: usize,
    pub n_samples: usize,
    pub n_grm_elements: usize,
    pub mean_kinship: f64,
    pub sd_kinship: f64,
    pub min_kinship: f64,
    pub max_kinship: f64,
}