#!usr/bin/R

# convert_hg19to38_repli_time.R
# make replication timing file for ACE-Seq usage by converting hg19 to hg38 coordinates
# usage: Rscript convert_hg19to38_repli_time.R <hg19_aceseq_rt_file> <hg19_to_hg38_chain> <outdir>

# load libraries
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

hg19_aceseq_rt_file <- args[1]
hg19_to_hg38_chain <- args[2]
outdir <- args[3]

# read in hg19 RT file ReplicationTime_10cellines_mean_10KB.Rda
load(hg19_aceseq_rt_file)
if (!exists("time10")) {
    stop("ReplicationTime_10cellines_mean_10KB.Rda was not found.")
}

# convert to GRanges
chroms <- as.character(time10$CHROM)
if (!any(grepl("chr", chroms[1:10]))) {
    chroms <- paste0("chr", chroms)
}

## convert 10kb index to start/end
## tenkb 2 = 10,001 - 20,000
starts <- (time10$tenkb - 1) * 10000 + 1
ends   <- time10$tenkb * 10000

gr_hg19 <- GRanges(
    seqnames = chroms,
    ranges   = IRanges(start = starts, end = ends),
    time     = time10$time
)

# import chain file
chain <- import.chain(hg19_to_hg38_chain)
lifted_list <- liftOver(gr_hg19, chain)
gr_hg38_raw <- unlist(lifted_list)

# tile the genome into 10kb windows
sizes_df <- data.frame(
    chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
            "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
            "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),
    len = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 
            159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 
            114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 
            58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
)
seq_info <- Seqinfo(seqnames = sizes_df$chr, seqlengths = sizes_df$len, genome = "hg38")
gr_hg38_grid <- tileGenome(seq_info, tilewidth = 10000, cut.last.tile.in.chrom = TRUE)

# map lifted data to the new grid
hits <- findOverlaps(gr_hg38_grid, gr_hg38_raw)

# For every new bin, take the value of the old bin that overlaps it.
# If multiple old bins overlap one new bin, we take the one with the *maximum* overlap width.

map_df <- data.frame(
    queryIdx = queryHits(hits), 
    subjectIdx = subjectHits(hits),
    overlap_width = width(pintersect(gr_hg38_grid[queryHits(hits)], gr_hg38_raw[subjectHits(hits)]))
)

# sort by overlap width so the best overlap comes first
map_df <- map_df %>% 
    arrange(desc(overlap_width)) %>%
    filter(!duplicated(queryIdx)) # Keep only the best match for each new bin

gr_hg38_grid$time <- NA
gr_hg38_grid$time[map_df$queryIdx] <- gr_hg38_raw$time[map_df$subjectIdx]

# remove bins that got no data (NAs)
gr_hg38_final <- gr_hg38_grid[!is.na(gr_hg38_grid$time)]

# convert back to the dataframe format ACEseq uses: CHROM, tenkb, time
# tenkb index = floor(end / 10000)
# ACEseq hg38 expects 'chr1', 'chr2' etc.
output_df <- data.frame(
    CHROM = as.character(seqnames(gr_hg38_final)),
    tenkb = ceiling(end(gr_hg38_final) / 10000),
    time  = gr_hg38_final$time
)

time10_hg38 <- output_df

# save
save(time10_hg38, file = file.path(outdir, "ReplicationTime_10cellines_mean_10KB_hg38.Rda"))
cat("Done! Saved as ReplicationTime_10cellines_mean_10KB_hg38.Rda\n")