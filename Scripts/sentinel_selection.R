library(data.table)

# config 

traits <- c("LD_ClUMPED_FINAL")

input_base <- "~/Desktop/GWAS_Files"
output_dir <- "~/Desktop/GWAS_Files"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

WIDTH       <- 1e6       # ±1 Mb
PTHRESH      <- 5e-9
chromosomes <- 1:22

# expected columns

clean_cols <- c(
  "SNP",
  "chromosome",
  "position",
  "p_value_lrt",
  "beta1",
  "se1"
)

# sentinel selection

selectSentinels <- function(dt, width_bp = 1e6) {
  dt <- copy(dt)
  setorder(dt, pvalue)
  
  sent <- dt[0]
  
  while (nrow(dt) > 0) {
    lead <- dt[1]
    
    sent <- rbind(sent, lead, use.names = TRUE, fill = TRUE)
    
    dt <- dt[
      chromosome != lead$chromosome |
        position < (lead$position - width_bp) |
        position > (lead$position + width_bp)
    ]
  }
  
  sent
}

# loop 

for (trait in traits) {
  
  infile <- file.path(
    input_base,
    sprintf("%s.clumps", trait)
  )
  
  if (!file.exists(infile)) {
    warning("File not found: ", infile)
    next
  }
  
  message("Processing trait: ", trait)
  message("Input file: ", infile)
  
  
  dt <- fread(
    infile,
    select = c("rs_number", "p-value", "beta", "se", "chr", "pos")
  )
  
  # rename columns as per file
  setnames(
    dt, 
    old = c("rs_number", "chr", "pos", "p-value", "beta", "se"), 
    new = c("SNP", "chromosome", "position", "p_value_lrt", "beta1", "se1")
  )
  
  
  dt[, chromosome := as.integer(chromosome)]
  dt[, position   := as.integer(position)]
  dt[, pvalue     := as.numeric(p_value_lrt)]
  
  # autosomes only
  dt <- dt[chromosome %in% chromosomes]
  
  # filter p-value threshold
  dt <- dt[p_value_lrt < PTHRESH]
  
  if (!nrow(dt)) {
    message("No variants passed p < ", PTHRESH, " for ", trait)
    next
  }
  
  # Sentinel selection across the whole trait file
  sentinels <- selectSentinels(dt, WIDTH)
  
  sentinels_clean <- sentinels[, ..clean_cols]
  
  setorder(sentinels_clean, chromosome, position)
  
  outfile <- file.path(
    output_dir,
    sprintf("%s_sentinels_clean.txt", trait)
  )
  
  fwrite(sentinels_clean, outfile, sep = "\t", quote = FALSE)
  
  message("Sentinel file written to: ", outfile)
}