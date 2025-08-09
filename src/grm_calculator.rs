
use crate::data_io::VntrData;
use crate::error::Result;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use ndarray::s;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::sync::{Arc, Mutex};
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
        
        // Step 1: Calculate allele frequencies for each variant (using a dedicated rayon pool)
        log::info!("Calculating allele frequencies for {} variants...", n_variants);
        let pool = ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .expect("Failed to build rayon thread pool");
        let allele_freqs: Vec<f64> = pool.install(|| {
            (0..n_variants)
                .into_par_iter()
                .map(|variant_idx| {
                    let mut sum = 0.0;
                    let mut valid_count = 0;
                    for sample_idx in 0..n_samples {
                        let dosage = dosage_matrix[[sample_idx, variant_idx]];
                        if dosage.is_finite() {
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
                .collect()
        });
        
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

        // Build standardized matrix Z with empirical per-variant mean/std
        // First compute empirical stats and filter low-variance columns
        let std_epsilon: f64 = 1e-8;
        log::info!("Standardization: empirical per-variant mean/std with variance filtering (std >= {}, count >= 2)", std_epsilon);
        pb.set_message("Computing per-variant mean/std");
        let variant_stats_all: Vec<(f64, f64, usize)> = valid_variants
            .par_iter()
            .map(|&variant_idx| {
                let mut sum = 0.0f64;
                let mut sumsq = 0.0f64;
                let mut count: usize = 0;
                for sample_idx in 0..n_samples {
                    let v = dosage_matrix[[sample_idx, variant_idx]];
                    if v.is_finite() {
                        sum += v;
                        sumsq += v * v;
                        count += 1;
                    }
                }
                if count >= 2 {
                    let mean = sum / count as f64;
                    let var = (sumsq - sum * sum / count as f64) / (count - 1) as f64;
                    (mean, var.sqrt(), count)
                } else if count == 1 {
                    // Single observation
                    let mean = sum;
                    (mean, 0.0, count)
                } else {
                    (0.0, 0.0, 0)
                }
            })
            .collect();

        // Filter effective variants by std and count
        let mut effective_indices: Vec<usize> = Vec::new();
        let mut variant_stats: Vec<(f64, f64)> = Vec::new();
        for (pos, &vidx) in valid_variants.iter().enumerate() {
            let (mean_k, std_k, cnt_k) = variant_stats_all[pos];
            if cnt_k >= 2 && std_k >= std_epsilon {
                effective_indices.push(vidx);
                variant_stats.push((mean_k, std_k));
            }
        }
        let m_effective = effective_indices.len();
        if m_effective == 0 {
            return Err(crate::error::DgrmError::InvalidFormat(
                "No informative variants after variance filtering".to_string()
            ));
        }
        log::info!("Effective variants after variance filter: {} (from {} after MAF)", m_effective, n_valid_variants);

        pb.set_message("Building standardized matrix Z");
        let mut z = Array2::<f64>::zeros((n_samples, m_effective));
        for sample_idx in 0..n_samples {
            let mut row = z.row_mut(sample_idx);
            for (new_idx, &variant_idx) in effective_indices.iter().enumerate() {
                let v = dosage_matrix[[sample_idx, variant_idx]];
                if v.is_finite() {
                    let (mean_k, std_k) = variant_stats[new_idx];
                    row[new_idx] = (v - mean_k) / std_k;
                } else {
                    row[new_idx] = 0.0;
                }
            }
        }

        // Compute GRM = Z * Z^T / (m-1) using upper-triangular blocked parallel matmul
        pb.set_message("Matrix multiply Z * Z^T (blocked upper-triangular)");
        let z_t = z.view().reversed_axes();
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

        use std::sync::Mutex;
        let grm_matrix_arc = Arc::new(Mutex::new(Array2::<f64>::zeros((n_samples, n_samples))));
        pool.install(|| {
            coords.into_par_iter().for_each(|(rs, re, cs, ce)| {
                let a = z.slice(s![rs..re, ..]);
                let b = z_t.slice(s![.., cs..ce]);
                let sub = a.dot(&b);
                let mut guard = grm_matrix_arc.lock().unwrap();
                guard.slice_mut(s![rs..re, cs..ce]).assign(&sub);
                if cs != rs {
                    guard.slice_mut(s![cs..ce, rs..re]).assign(&sub.t());
                }
            });
        });
        let mut grm_matrix = Arc::try_unwrap(grm_matrix_arc)
            .expect("Multiple references to grm_matrix remain")
            .into_inner()
            .expect("Mutex poisoned");
        let scale = 1.0 / (m_effective.max(1)) as f64;
        grm_matrix.mapv_inplace(|x| x * scale);

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
            n_variants: m_effective, // Use effective variant count
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
        let grm_matrix = Arc::new(Mutex::new(Array2::<f64>::zeros((n_samples, n_samples))));
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
        
        // Progress tracking
        let pb = ProgressBar::new(total_block_pairs as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} blocks ({percent}%) {msg}")
            .unwrap()
            .progress_chars("#>-"));
        pb.set_message("Processing blocks");

        // Use empirical per-variant standardization; no per-sample stats needed
        log::info!("Standardization: empirical per-variant mean/std");

        // Anti-diagonal waves to ensure disjoint writes; parallelize within each wave
        for d in 0..n_blocks {
            let pairs: Vec<(usize, usize)> = (0..(n_blocks - d)).map(|i| (i, i + d)).collect();
            pairs.into_par_iter().for_each(|(block_i, block_j)| {
                let start_i = block_i * block_size;
                let end_i = (start_i + block_size).min(n_samples);
                let start_j = block_j * block_size;
                let end_j = (start_j + block_size).min(n_samples);

                // Build Z blocks for I and J using per-variant standardization
                let rows_i = end_i - start_i;
                let rows_j = end_j - start_j;
                // Pre-compute empirical mean/std per variant once and filter low-variance
                let std_epsilon: f64 = 1e-8;
                let stats_all: Vec<(f64, f64, usize)> = valid_variants
                    .iter()
                    .map(|&variant_idx| {
                        let mut sum = 0.0f64;
                        let mut sumsq = 0.0f64;
                        let mut count: usize = 0;
                        for i in 0..n_samples {
                            let v = dosage_matrix[[i, variant_idx]];
                            if v.is_finite() {
                                sum += v;
                                sumsq += v * v;
                                count += 1;
                            }
                        }
                        if count >= 2 {
                            let mean = sum / count as f64;
                            let var = (sumsq - sum * sum / count as f64) / (count - 1) as f64;
                            (mean, var.sqrt(), count)
                        } else if count == 1 {
                            let mean = sum;
                            (mean, 0.0, count)
                        } else {
                            (0.0, 0.0, 0)
                        }
                    })
                    .collect();

                let mut eff_indices: Vec<usize> = Vec::new();
                let mut stats: Vec<(f64, f64)> = Vec::new();
                for (pos, &vidx) in valid_variants.iter().enumerate() {
                    let (mean_k, std_k, cnt_k) = stats_all[pos];
                    if cnt_k >= 2 && std_k >= std_epsilon {
                        eff_indices.push(vidx);
                        stats.push((mean_k, std_k));
                    }
                }

                let m = eff_indices.len();
                let mut z_i = Array2::<f64>::zeros((rows_i, m));
                let mut z_j = Array2::<f64>::zeros((rows_j, m));

                for (k, &variant_idx) in eff_indices.iter().enumerate() {
                    let (mean_k, std_k) = stats[k];
                    for (local_r, i) in (start_i..end_i).enumerate() {
                        let x = dosage_matrix[[i, variant_idx]];
                        z_i[[local_r, k]] = if x.is_finite() { (x - mean_k) / std_k } else { 0.0 };
                    }
                    for (local_r, j) in (start_j..end_j).enumerate() {
                        let x = dosage_matrix[[j, variant_idx]];
                        z_j[[local_r, k]] = if x.is_finite() { (x - mean_k) / std_k } else { 0.0 };
                    }
                }

                // Multiply: sub = Z_I * Z_J^T
                let sub = z_i.dot(&z_j.t());

                // Write upper and mirror to lower if needed
                let mut guard = grm_matrix.lock().unwrap();
                guard
                    .slice_mut(s![start_i..end_i, start_j..end_j])
                    .assign(&sub);
                if block_i != block_j {
                    guard
                        .slice_mut(s![start_j..end_j, start_i..end_i])
                        .assign(&sub.t());
                }

                pb.inc(1);
            });
        }

        pb.finish_with_message("Block-wise calculation completed");
        
        let calc_time = start_time.elapsed();
        log::info!("Block-wise GRM calculation completed in {:.2} seconds", calc_time.as_secs_f64());

        // Scale by 1/m (no diagonal normalization)
        let mut grm_matrix = Arc::try_unwrap(grm_matrix)
            .expect("Multiple references to grm_matrix remain")
            .into_inner()
            .expect("Mutex poisoned");
        let scale = 1.0 / (valid_variants.len().max(1)) as f64;
        grm_matrix.mapv_inplace(|x| x * scale);
        
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
            n_variants: valid_variants.len(),
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