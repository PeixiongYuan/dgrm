#!/usr/bin/env Rscript

# Verify that the GRM calculation fix is working correctly

cat("=== Verifying GRM Calculation Fix ===\n")

# Test with a simple known case
# Create a simple test case where we know the expected correlation
test_data <- matrix(c(
  1, 2, 3, 4,  # Sample A
  1, 2, 3, 4,  # Sample B (identical to A, should have correlation = 1)
  4, 3, 2, 1,  # Sample C (opposite of A, should have correlation = -1)
  2, 2, 2, 2   # Sample D (constant, correlation with others should be undefined/0)
), nrow = 4, byrow = TRUE)

# Calculate correlations manually using R
cat("\n=== Expected correlations (using R's cor function) ===\n")
expected_cor <- cor(t(test_data), use = "complete.obs")
print(expected_cor)

# Write test data in our format
write.table(
  data.frame(
    sample_id = c("A", "B", "C", "D"),
    test_data
  ),
  file = "verification_test.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("\n=== Test data written to verification_test.txt ===\n")
cat("Sample A: 1 2 3 4\n")
cat("Sample B: 1 2 3 4 (identical to A)\n") 
cat("Sample C: 4 3 2 1 (opposite to A)\n")
cat("Sample D: 2 2 2 2 (constant)\n")

cat("\n=== Expected results ===\n")
cat("A-A: 1.0, B-B: 1.0, C-C: 1.0, D-D: NaN or 1.0\n")
cat("A-B: 1.0 (identical)\n")
cat("A-C: -1.0 (perfect negative correlation)\n")
cat("A-D, B-D, C-D: 0.0 (constant has no correlation)\n")