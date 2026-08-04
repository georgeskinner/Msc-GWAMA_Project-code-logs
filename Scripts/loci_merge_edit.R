sink(stderr())
library(data.table)
library(dplyr)
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))


traits <- c("LD_CLUMPED_FINAL") 

input_base  <- "~/Desktop/GWAS_Files"
output_dir  <- "~/Desktop/GWAS_Files"

# parameters
width <- 1e6   # +/- 1Mb window
gap   <- 0     # merge max gap

# loop
for (trait in traits) {
  
  input_file  <- file.path(input_base, sprintf("%s_sentinels_clean.txt", trait))
  output_file <- file.path(output_dir, sprintf("%s_sentinel_merge_loci.txt", trait))
  
  if (!file.exists(input_file)) {
    warning("Sentinel file not found: ", input_file)
    next
  }
  
  message("Processing trait merging for: ", trait)
  
  # load sentinels
  sentinels <- fread(input_file) %>%
    arrange(chromosome, position) %>%
    mutate(
      start = position - width,
      end   = position + width
    )
  
  if (nrow(sentinels) == 0) {
    message("No sentinels found in file: ", input_file)
    next
  }
  
  
  gr_sentinels <- makeGRangesFromDataFrame(
    sentinels,
    seqnames.field = "chromosome",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE
  )
  
  
  gr_loci <- GenomicRanges::reduce(gr_sentinels, min.gapwidth = gap)
  
  
  loci_df <- as_tibble(gr_loci) %>%
    mutate(locus = paste0("locus", row_number())) %>%
    dplyr::select(locus, chromosome = seqnames, locus_start = start, locus_end = end, locus_size = width)
  
  
  overlaps <- findOverlaps(gr_sentinels, gr_loci)
  sentinels$locus <- loci_df$locus[subjectHits(overlaps)]
  
  
  loci_summary <- sentinels %>%
    group_by(locus) %>%
    arrange(pvalue) %>% 
    summarise(
      lead_SNP    = dplyr::first(SNP),
      lead_p      = dplyr::first(pvalue),
      .groups     = "drop"
    )
  
  # join genomic coordinates back to the summary
  final_loci <- inner_join(loci_df, loci_summary, by = "locus") %>%
    arrange(as.integer(chromosome), locus_start)
  
  # save
  write.table(final_loci, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("✅ Locus merge table written to: ", output_file)
}