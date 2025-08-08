
use crate::data_io::VntrData;
use crate::error::Result;
use crate::grm_calculator::{GrmMatrix, GrmStats};
use byteorder::{LittleEndian, WriteBytesExt};
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
    pub fn write_grm(&self, grm_matrix: &GrmMatrix, _vntr_data: &VntrData) -> Result<()> {
        log::info!("Writing GCTA GRM files with prefix: {}", self.output_prefix);

        // Calculate statistics
        let stats = grm_matrix.calculate_stats();

        // Write binary files
        self.write_grm_n_bin(grm_matrix)?;
        self.write_grm_bin(grm_matrix)?;
        
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

    /// Write number of variants per GRM element (4-byte integer for each upper-triangular element)
    /// GCTA expects one 4-byte integer per element in the same order as .grm.bin
    fn write_grm_n_bin(&self, grm_matrix: &GrmMatrix) -> Result<()> {
        let file_path = format!("{}.grm.N.bin", self.output_prefix);
        let mut file = File::create(&file_path)?;

        let n_samples = grm_matrix.n_samples();
        let n_variants = grm_matrix.n_variants();
        let num_elements = (n_samples * (n_samples + 1)) / 2; // upper triangular with diagonal

        for _ in 0..num_elements {
            file.write_i32::<LittleEndian>(n_variants as i32)?;
        }

        log::debug!(
            "Written N per-element ({} ints, value = {}) to {}",
            num_elements,
            n_variants,
            file_path
        );
        Ok(())
    }

    /// Write GRM values as 4-byte floats (upper triangular including diagonal)
    fn write_grm_bin(&self, grm_matrix: &GrmMatrix) -> Result<()> {
        let file_path = format!("{}.grm.bin", self.output_prefix);
        let file = File::create(&file_path)?;
        let mut writer = BufWriter::new(file);
        
        let values = grm_matrix.upper_triangular_values();
        
        for &value in &values {
            writer.write_f32::<LittleEndian>(value as f32)?;
        }
        
        writer.flush()?;
        
        log::debug!("Written {} GRM values to {}", values.len(), file_path);
        Ok(())
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