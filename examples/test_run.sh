#!/bin/bash

# Test script for DGRM VNTR GRM Calculator

set -e

echo "=== DGRM Test Run ==="

# Check if binary exists
if [ ! -f "../target/release/dgrm" ]; then
    echo "Binary not found. Building first..."
    cd ..
    cargo build --release
    cd examples
fi

# Create test directory
mkdir -p test_output
cd test_output

# Run with sample data
echo "Running DGRM with sample data..."
../../target/release/dgrm \
    --input ../sample_data.txt \
    --output sample_grm \
    --threads 4 \
    --verbose

# Check output files
echo ""
echo "=== Output Files ==="
ls -la sample_grm.*

echo ""
echo "=== File Contents ==="

echo "Sample IDs (first 5 lines):"
head -5 sample_grm.grm.id

echo ""
echo "Summary statistics:"
cat sample_grm.grm.summary

echo ""
echo "File sizes:"
du -h sample_grm.*

echo ""
echo "=== Test Completed Successfully ==="
echo "GCTA format files are ready for use with GCTA software."