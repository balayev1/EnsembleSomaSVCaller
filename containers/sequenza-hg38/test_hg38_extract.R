#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: test_hg38_extract.R <merged.seqz.gz>", call. = FALSE)
}

seqz_file <- args[[1]]

library(sequenza)
library(copynumber)

Sys.setenv(VROOM_CONNECTION_SIZE = Sys.getenv("VROOM_CONNECTION_SIZE", unset = "33554432"))

cat("sequenza:", as.character(packageVersion("sequenza")), "\n")
cat("copynumber:", as.character(packageVersion("copynumber")), "\n")
cat("seqz:", seqz_file, "\n")

result <- sequenza.extract(
    file = seqz_file,
    assembly = "hg38",
    verbose = FALSE,
    parallel = 1
)

cat("extract fields:", paste(names(result), collapse = ", "), "\n")
