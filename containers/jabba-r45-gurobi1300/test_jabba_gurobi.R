#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(JaBbA)
  library(slam)
  library(gurobi)
  library(GenomicRanges)
  library(IRanges)
  library(gUtils)
})

cat("JaBbA:", as.character(packageVersion("JaBbA")), "\n")
cat("slam:", as.character(packageVersion("slam")), "\n")
cat("gurobi:", as.character(packageVersion("gurobi")), "\n")
cat("GUROBI_HOME:", Sys.getenv("GUROBI_HOME", unset = "<unset>"), "\n")
cat("GRB_LICENSE_FILE:", Sys.getenv("GRB_LICENSE_FILE", unset = "<unset>"), "\n")

gr <- GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(c(1L, 5L), width = 1L))
cat("seqnames(GRanges):", paste(as.character(seqnames(gr)), collapse = ","), "\n")

grl <- split(gr, c("pair", "pair"))
grl_nochr <- gUtils::gr.nochr(grl)
cat("seqlevels(gr.nochr(GRangesList)):", paste(seqlevels(grl_nochr), collapse = ","), "\n")
