
use crate::data_io::VntrData;
use crate::error::Result;
use crate::grm_calculator::{GrmMatrix, GrmStats};
use bytemuck::cast_slice;
use std::fs::File;
use std::io::{BufWriter, Write};


pub struct GctaWriter {
    output_prefix: String,
}

impl GctaWriter {
    pub fn new(output_prefix: &str) -> Self {
        Self {
            output_prefix: output_prefix.to_string(),
        }
    }

    /// Write GRM in GCTA format
    pub fn write_grm(&self, grm_matrix: &GrmMatrix, vntr_data: &VntrData) -> Result<()> {
        log::info!("Writing GCTA GRM files with prefix: {}", self.output_prefix);
        
        // Write binary files (streaming I/O), compute stats while writing GRM values
        self.write_grm_n_bin(grm_matrix, vntr_data)?;
        let stats = self.write_grm_bin(grm_matrix)?;
        
        // Write text files
        self.write_grm_id(grm_matrix.sample_ids())?;
        self.write_grm_summary(&stats)?;

        log::info!("GCTA GRM files created successfully:");
        log::info!("  1. {}.grm.N.bin - Number of variants (binary)", self.output_prefix);
        log::info!("  2. {}.grm.bin - GRM values (binary)", self.output_prefix);
        log::info!("  3. {}.grm.id - Sample IDs (text)", self.output_prefix);
        log::info!("  4. {}.grm.summary - Summary statistics (text)", self.output_prefix);

        self.print_summary(&stats);

        Ok(())
    }

    /// Write number of variants per GRM element (4-byte FLOAT per element), row-major lower triangular
    /// Stream in chunks to avoid allocating O(n^2) buffer
    fn write_grm_n_bin(&self, grm_matrix: &GrmMatrix, _vntr_data: &VntrData) -> Result<()> {
        let file_path = format!("{}.grm.N.bin", self.output_prefix);
        let file = File::create(&file_path)?;
        let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, file);

        let n = grm_matrix.n_samples();
        let m = grm_matrix.n_variants() as f32;
        let len = n * (n + 1) / 2;
        // Stream write constant value m
        let chunk_elems: usize = 4_194_304; // ~16MB per chunk
        let chunk: Vec<f32> = vec![m; chunk_elems];
        let mut remaining = len;
        while remaining >= chunk_elems {
            writer.write_all(cast_slice(&chunk))?;
            remaining -= chunk_elems;
        }
        if remaining > 0 {
            writer.write_all(cast_slice(&chunk[..remaining]))?;
        }

        writer.flush()?;
        log::debug!("Written per-pair N (float) in row-major lower-triangular order to {}", file_path);
        Ok(())
    }

    /// Write GRM values as 4-byte floats, row-major lower-triangular (including diagonal)
    /// Stream per-row to avoid allocating O(n^2) buffer, and compute stats on the fly
    fn write_grm_bin(&self, grm_matrix: &GrmMatrix) -> Result<GrmStats> {
        let file_path = format!("{}.grm.bin", self.output_prefix);
        let file = File::create(&file_path)?;
        let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, file);

        let n = grm_matrix.n_samples();
        let mut count: usize = 0;
        let mut mean: f64 = 0.0;  // Welford online mean
        let mut m2: f64 = 0.0;    // Welford accumulated variance numerator
        let mut min_val: f64 = f64::INFINITY;
        let mut max_val: f64 = f64::NEG_INFINITY;

        // Row-major lower triangular streaming
        for j in 0..n {
            let mut row_buf: Vec<f32> = Vec::with_capacity(j + 1);
            for k in 0..=j {
                let v = grm_matrix.value(j, k);
                // stats
                count += 1;
                let delta = v - mean;
                mean += delta / (count as f64);
                m2 += delta * (v - mean);
                if v < min_val { min_val = v; }
                if v > max_val { max_val = v; }
                // buffer
                row_buf.push(v as f32);
            }
            writer.write_all(cast_slice(&row_buf))?;
        }

        writer.flush()?;
        log::debug!("Written GRM values (float) in row-major lower-triangular order to {}", file_path);

        let n_elements = count;
        let variance = if n_elements > 1 { m2 / (n_elements as f64) } else { 0.0 };
        let stats = GrmStats {
            n_variants: grm_matrix.n_variants(),
            n_samples: n,
            n_grm_elements: n_elements,
            mean_kinship: mean,
            sd_kinship: variance.sqrt(),
            min_kinship: min_val,
            max_kinship: max_val,
        };
        Ok(stats)
    }

    /// Write sample IDs in GCTA format (FID\tIID)
    fn write_grm_id(&self, sample_ids: &[String]) -> Result<()> {
        let file_path = format!("{}.grm.id", self.output_prefix);
        let file = File::create(&file_path)?;
        let mut writer = BufWriter::new(file);
        
        for sample_id in sample_ids {
            writeln!(writer, "{}\t{}", sample_id, sample_id)?;
        }
        
        writer.flush()?;
        
        log::debug!("Written {} sample IDs to {}", sample_ids.len(), file_path);
        Ok(())
    }

    /// Write summary statistics
    fn write_grm_summary(&self, stats: &GrmStats) -> Result<()> {
        let file_path = format!("{}.grm.summary", self.output_prefix);
        let file = File::create(&file_path)?;
        let mut writer = BufWriter::new(file);
        
        writeln!(writer, "Metric\tValue")?;
        writeln!(writer, "N_variants\t{}", stats.n_variants)?;
        writeln!(writer, "N_samples\t{}", stats.n_samples)?;
        writeln!(writer, "N_GRM_elements\t{}", stats.n_grm_elements)?;
        writeln!(writer, "Mean_kinship\t{:.6}", stats.mean_kinship)?;
        writeln!(writer, "SD_kinship\t{:.6}", stats.sd_kinship)?;
        writeln!(writer, "Min_kinship\t{:.6}", stats.min_kinship)?;
        writeln!(writer, "Max_kinship\t{:.6}", stats.max_kinship)?;
        
        writer.flush()?;
        
        log::debug!("Written summary statistics to {}", file_path);
        Ok(())
    }

    /// Print summary statistics to console
    fn print_summary(&self, stats: &GrmStats) {
        log::info!("Summary statistics:");
        log::info!("  N_variants: {}", stats.n_variants);
        log::info!("  N_samples: {}", stats.n_samples);
        log::info!("  N_GRM_elements: {}", stats.n_grm_elements);
        log::info!("  Mean_kinship: {:.6}", stats.mean_kinship);
        log::info!("  SD_kinship: {:.6}", stats.sd_kinship);
        log::info!("  Min_kinship: {:.6}", stats.min_kinship);
        log::info!("  Max_kinship: {:.6}", stats.max_kinship);
        log::info!("Files are ready for use with GCTA software.");
    }


}