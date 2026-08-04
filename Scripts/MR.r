library(TwoSampleMR)
library(data.table)

data_folder <- "/Users/George/Desktop/GWAS_Files/FinnGen_Downloads"
exposure_dat <- readRDS("FINAL_EXPOSURE_PREPARED.rds")


outcome_files <- list.files(path = data_folder, pattern = "finngen_R10_.*\\.gz$", full.names = TRUE)

results_list <- list()

# Loop through each FinnGen file
for (file in outcome_files) {
  message("\n--- Processing: ", basename(file), " ---")
  
  tryCatch({
    # Read outcome data
    outcome_dat <- read_outcome_data(
      filename = file,
      snps = exposure_dat$SNP,
      sep = "\t",
      snp_col = "rsids",
      beta_col = "beta",
      se_col = "sebeta",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      pval_col = "pval",
      eaf_col = "af_alt" # Critical for harmonizing palindromic SNPs
    )
    
    # Check if outcome_dat is empty before proceeding
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      message("Skipping: No matching SNPs found in ", basename(file))
      next
    }
    
    # Harmonize (action = 2 is default, handles palindromic SNPs safely)
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    
    # Check if any SNPs survived harmonization
    if (nrow(dat) > 0) {
      res <- mr(dat)
      
      # Ensure mr() actually produced a result before trying to add the filename
      if (!is.null(res) && nrow(res) > 0) {
        res$outcome_file <- basename(file)
        results_list[[basename(file)]] <- res
        message("Successfully processed: ", basename(file))
      } else {
        message("Not enough valid SNPs to perform MR for ", basename(file))
      }
      
    } else {
      message("No SNPs remained after harmonization for ", basename(file))
    }
    
  }, error = function(e) {
    message("Error processing ", basename(file), ": ", e$message)
  })
}

# save results

if (length(results_list) > 0) {
  final_results <- do.call(rbind, results_list)
  write.csv(final_results, "MR_Results.csv", row.names = FALSE)
  message("\nDone! Results saved to 'MR_Results.csv'.")
} else {
  message("\nNo results to save. All files failed or had no matching SNPs.")
}


