
library(GenomicSEM)

# config 

gwama_file       <- "GWAMA_FINAL_FIXED_COPY.out"
gwama_name       <- "GWAMA_Primary_Trait"
gwama_N          <- 2310139  
gwama_samp_prev  <- NA     
gwama_pop_prev   <- NA     

finngen_traits <- list(
  list(file = "FinnGen_Downloads/finngen_R10_CD2_BENIGN_LIVER.gz",            name = "Benign_Liver",                 samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_K11_LIVER.gz",                  name = "Liver_Disease",                samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_K11_IBD_STRICT.gz",             name = "IBD_Strict",                   samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_RHEU_ARTHRITIS_OTH.gz",         name = "Rheumatoid_Arthritis_Other",   samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_N14_CHRONKIDNEYDIS.gz",         name = "Chronic_Kidney_Disease",       samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_I9_HEARTFAIL_AND_OVERWEIGHT.gz", name = "Heart_Failure_And_Overweight", samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_COPD_OPPORTUNIST_INFECTIONS.gz", name = "COPD_Opportunist_Infections",  samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_COPD_LATER.gz",                 name = "COPD_Later",                   samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_COPD_INSUFFICIENCY.gz",          name = "COPD_Insufficiency",           samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_COPD_EARLY.gz",                  name = "COPD_Early",                   samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_M13_RHEUMA.gz",                 name = "Rheuma",                       samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_G6_ALZHEIMER.gz",               name = "Alzheimers",                   samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_J10_COPD.gz",                   name = "COPD_General",                 samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_I9_STR.gz",                     name = "Stroke",                       samp_prev = NA, pop_prev = NA),
  list(file = "FinnGen_Downloads/finngen_R10_I9_HEARTFAIL.gz",               name = "Heart_Failure",                samp_prev = NA, pop_prev = NA)
)


# munge for primary trait

cat("\n--- Munging Primary GWAMA Trait ---\n")
gwama_munged <- paste0(gwama_name, ".sumstats.gz")

if (!file.exists(gwama_munged)) {
  munge(
    files = gwama_file,
    hm3 = "w_hm3.snplist",
    trait.names = gwama_name,
    N = gwama_N,
    info.filter = 0.9,
    maf.filter = 0.01,
    column.names = list(
      SNP = "rs_number",
      A1 = "reference_allele",
      A2 = "other_allele",
      effect = "beta",
      P = "p",
      MAF = "eaf"
    )
  )
} else {
  cat("Munged file already exists. Skipping munge for Primary Trait.\n")
}

# finngen loop

results_list <- list()

for (i in seq_along(finngen_traits)) {
  fg <- finngen_traits[[i]]
  cat(sprintf("\n=== Processing [%d/%d]: %s ===\n", i, length(finngen_traits), fg$name))
  
  fg_munged <- paste0(fg$name, ".sumstats.gz")
  tmp_clean <- paste0("temp_clean_", fg$name, ".txt")
  
  if (!file.exists(fg_munged)) {
    tryCatch({
      cat("Safely loading and standardizing FinnGen file for", fg$name, "...\n")
      
     
      dt <- data.table::fread(fg$file, header = FALSE, skip = 1, tmpdir = getwd())
      
      
      dt_clean <- dt[, .(
        SNP    = V5,
        A1     = V4,
        A2     = V3,
        effect = V9,
        P      = V7,
        MAF    = V11,
        N      = 412181  
      )]
      
      
      data.table::fwrite(dt_clean, tmp_clean, sep = "\t")
      rm(dt, dt_clean); gc() 
      
      cat("Munging standardized file...\n")
      munge(
        files = tmp_clean,
        hm3 = "w_hm3.snplist",
        trait.names = fg$name,
        info.filter = 0.9,
        maf.filter = 0.01
      )
      
      # Delete temporary clean file
      if (file.exists(tmp_clean)) file.remove(tmp_clean)
      
    }, error = function(e) {
      cat("Error preparing/munging file:", fg$name, "\n", conditionMessage(e), "\n")
      if (file.exists(tmp_clean)) file.remove(tmp_clean)
    })
  }
  
  # run bivariate LDSC
  
  if (file.exists(fg_munged)) {
    tryCatch({
      ldsc_out <- ldsc(
        traits = c(gwama_munged, fg_munged),
        sample.prev = c(gwama_samp_prev, fg$samp_prev),
        population.prev = c(gwama_pop_prev, fg$pop_prev),
        ld = "eur_w_ld_chr/",
        wld = "eur_w_ld_chr/"
      )
      
      cov_g  <- ldsc_out$S[1, 2]
      var_g1 <- ldsc_out$S[1, 1]
      var_g2 <- ldsc_out$S[2, 2]
      
      rg    <- cov_g / sqrt(var_g1 * var_g2)
      rg_se <- sqrt(ldsc_out$V[2, 2])
      z_val <- rg / rg_se
      p_val <- 2 * pnorm(-abs(z_val))
      
      results_list[[i]] <- data.frame(
        Primary_Trait = gwama_name,
        FinnGen_Trait = fg$name,
        rg            = round(rg, 4),
        SE            = round(rg_se, 4),
        Z             = round(z_val, 3),
        P             = p_val,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      cat("Error running LDSC for:", fg$name, "\n", conditionMessage(e), "\n")
    })
  }
}


# p-value and save results

if (length(results_list) > 0) {
  final_results <- do.call(rbind, results_list)
  final_results$FDR_P <- p.adjust(final_results$P, method = "BH")
  
  cat("\n========================================================================\n")
  cat("                       BIVARIATE LDSC RESULTS                           \n")
  cat("========================================================================\n")
  print(final_results)
  
  write.csv(final_results, "LDSC_FinnGen_Correlations_Output.csv", row.names = FALSE)
  cat("\nResults successfully saved to 'LDSC_FinnGen_Correlations_Output.csv'\n")
} else {
  cat("\nNo results were generated. Please check for errors above.\n")
}