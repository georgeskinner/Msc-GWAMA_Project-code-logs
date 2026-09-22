library(data.table)
library(RMySQL)

gwama_res <- fread("GWAMA_FINAL_RANDOM.out", sep = "\t")


gwama_res[, chr := as.character(NA)]
gwama_res[, pos := as.integer(NA)]

setkey(gwama_res, rs_number)

# Extract unique RS IDs
rs_ids <- unique(gwama_res$rs_number)
rs_ids <- rs_ids[grepl("^rs[0-9]+$", rs_ids)]

cat(paste("Total unique RS IDs to map:", length(rs_ids), "\n"))

#Connect to UCSC
cat("Connecting to UCSC hg38 database stream...\n")
ucsc_db <- dbConnect(
  MySQL(),
  user = "genome",
  dbname = "hg38",
  host = "genome-mysql.soe.ucsc.edu",
  client.flag = 0 
)


chunk_size <- 50000
chunks <- split(rs_ids, ceiling(seq_along(rs_ids) / chunk_size))

for (i in seq_along(chunks)) {
  if (i %% 20 == 0) cat(sprintf("  Progress: Chunk %d of %d...\n", i, length(chunks)))
  
  
  id_list <- paste0("'", chunks[[i]], "'", collapse = ",")
  
  query <- paste0(
    "SELECT name AS rs_number, chrom AS CHR, chromEnd AS POS ",
    "FROM snp151 ",
    "WHERE name IN (", id_list, ")"
  )
  
 
  res <- dbGetQuery(ucsc_db, query)
  
  if (nrow(res) > 0) {
    tmp_dt <- as.data.table(res)
    

    tmp_dt[, CHR := gsub("^chr", "", CHR)]
    
    
    gwama_res[tmp_dt, on = .(rs_number), `:=`(chr = i.CHR, pos = as.integer(i.POS))]
  }
  
 
  if (i %% 50 == 0) gc()
}

dbDisconnect(ucsc_db)


cols <- colnames(gwama_res)
new_order <- c("rs_number", "chr", "pos", cols[!cols %in% c("rs_number", "chr", "pos")])
setcolorder(gwama_res, new_order)


fwrite(gwama_res, "GWAMA_FINAL_RANDOM_hg38.txt", sep = "\t", quote = FALSE, row.names = FALSE)