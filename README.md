## DGRM — Parallel GRM for VNTR

DGRM computes Genetic Relationship Matrices (GRM) from VNTR scaled dosage data. It is optimized for large cohorts, multi-threaded execution, and produces GCTA-compatible outputs.

## Highlights
- Multi-threaded core with automatic block-wise computation for large N
- Empirical per-variant standardization with low-variance filtering
- GCTA-compatible outputs: `.grm.bin`, `.grm.id`, `.grm.summary`, `.grm.N.bin`


## Build
```bash
cd dgrm
cargo build --release
./target/release/dgrm --help
```

## Usage
```bash
./target/release/dgrm \
  --input <FILE> \
  --output <PREFIX> \
  --threads <N> \
  [--keep <KEEP_FILE>]
```

Examples:
```bash
# Basic run
./target/release/dgrm --input chr1_vntr_scaled_dos.txt --output chr1_grm --threads 8

# Keep a subset of samples (FID IID or single-column file)
./target/release/dgrm --input data.txt --output out --threads 8 --keep keep.txt
```

## Input format
Tab-delimited, no header:
```
sample1	0.23	1.45	0.78	...
sample2	1.12	0.89	1.23	...
sample3	0.45	2.01	0.67	...
```
- Column 1: sample ID
- Columns 2..M: VNTR scaled dosage per variant (not necessarily 0/1/2)
- Rows are samples

## Method (concise)
- Empirical standardization per variant: for variant k, compute mean μ_k and sample std σ_k over finite values; set z[i,k] = (x[i,k] − μ_k)/σ_k; missing values become 0.
- Filter variants with low information: keep only σ_k ≥ 1e−8 and count ≥ 2; let m be the number of kept variants.
- Compute GRM: G = (Z · Zᵀ) / m. For very large N, G is built block-wise.

## Outputs (GCTA-compatible)
- `<prefix>.grm.bin`: GRM values as float32 in row-major lower-triangular order (including diagonal)
- `<prefix>.grm.N.bin`: number of variants used (float32, constant m for all pairs)
- `<prefix>.grm.id`: sample IDs (two columns: FID, IID)
- `<prefix>.grm.summary`: summary statistics

Note: This tool writes `.grm.bin` in the standard row-major lower-triangular order commonly accepted by GCTA utilities.

## Performance & memory
- Full GRM in memory uses ~ N² × 8 bytes (f64). The tool automatically switches to a block-wise algorithm for large N.
- Set threads via `--threads`. Throughput typically scales near physical core count.
- For extra speed, set `RUST_LOG=info` and consider tuning `DGRM_MM_CHUNK` (e.g., 256/512/1024).

## Logging
```bash
RUST_LOG=info ./target/release/dgrm --input data.txt --output result --threads 8 | cat
RUST_LOG=debug ./target/release/dgrm --input data.txt --output result --threads 8 | cat
```

## License
MIT