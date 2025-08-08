### DGRM — Parallel VNTR GRM Calculator (Rust)

DGRM is a high-performance Genetic Relationship Matrix (GRM) calculator for VNTR (Variable Number Tandem Repeat) scaled dosage data, implemented in Rust. It produces GCTA-compatible outputs and is optimized for multi-threaded execution and large cohorts.

### Features
- High-throughput, multi-threaded computation via Rayon
- Memory-aware block-wise algorithm for large N (auto-enabled on big inputs)
- GCTA-compatible outputs (.grm.bin, .grm.id, .grm.summary)
- Optional MAF filtering (min and/or max)
- Informative progress logging and summary statistics

### Build
```bash
cd dgrm
cargo build --release
./target/release/dgrm --help
```

### Usage
```bash
dgrm \
  --input <FILE> \
  --output <PREFIX> \
  --threads <N> \
  [--maf <FREQ>] [--max-maf <FREQ>]
```

Examples:
```bash
# Basic run
dgrm --input chr1_vntr_scaled_dos.txt --output chr1_grm --threads 8

# With MAF filtering
dgrm --input data.txt --output out --threads 16 --maf 0.01 --max-maf 0.99
```

### Input Format
Tab-delimited text without header:
```
sample1    0.23    1.45    0.78   ...
sample2    1.12    0.89    1.23   ...
sample3    0.45    2.01    0.67   ...
```
- First column: sample ID
- Remaining columns: scaled dosage per variant
- Rows = samples, columns = variants

### Output Files (GCTA-compatible)
- `<prefix>.grm.bin`: GRM values as float32 in column-major upper-triangular order (including diagonal)
- `<prefix>.grm.id`: sample IDs (two identical columns: FID, IID)
- `<prefix>.grm.summary`: summary statistics (text)
- `<prefix>.grm.N.bin`: number of variants used (4-byte integer)

Compatibility note:
- Some GCTA workflows expect `<prefix>.grm.N.bin` to contain per-pair counts in the same order as `<prefix>.grm.bin`. This tool currently writes a single total-variant count. The `.grm.bin` ordering is GCTA-conformant.

### Performance & Memory
- Full in-memory GRM requires ~N^2 × 8 bytes (float64). For very large N, the tool automatically switches to a memory-efficient block-wise correlation algorithm and writes into the final matrix incrementally.
- Ordering of values in `.grm.bin` matches GCTA: column-major, upper triangle (including diagonal).
- Set threads with `--threads`. Typical best throughput is near the number of physical cores.

### Logging
Use environment variable to increase verbosity:
```bash
RUST_LOG=info ./target/release/dgrm --input data.txt --output result --threads 8 | cat
RUST_LOG=debug ./target/release/dgrm --input data.txt --output result --threads 8 | cat
```

### License
MIT

### Version
0.1.0