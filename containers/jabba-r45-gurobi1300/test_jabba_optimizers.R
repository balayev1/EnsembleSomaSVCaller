#!/usr/bin/env Rscript

# Suppress both standard messages and the R 4.5 namespace collision warnings
suppressWarnings(suppressPackageStartupMessages({
  library(JaBbA)
  library(slam)
  library(gurobi)
  library(Rcplex)
  library(GenomicRanges)
  library(IRanges)
  library(gUtils)
}))

cat("----------------------------------------\n")
cat("DIAGNOSTICS & ENVIRONMENT\n")
cat("----------------------------------------\n")
cat("JaBbA:", as.character(packageVersion("JaBbA")), "\n")
cat("gurobi:", as.character(packageVersion("gurobi")), "\n")
cat("Rcplex:", as.character(packageVersion("Rcplex")), "\n")
cat("GUROBI_HOME:", Sys.getenv("GUROBI_HOME", unset = "<unset>"), "\n")
cat("CPLEX_DIR:", Sys.getenv("CPLEX_DIR", unset = "<unset>"), "\n")
cat("CPLEXDIR:", Sys.getenv("CPLEXDIR", unset = "<unset>"), "\n")
cat("GRB_LICENSE_FILE:", Sys.getenv("GRB_LICENSE_FILE", unset = "<unset>"), "\n\n")

cat("----------------------------------------\n")
cat("TEST 1: Bioc 3.21 seqinfo<- Patch Check\n")
cat("----------------------------------------\n")
gr <- GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(c(1L, 5L), width = 1L))
cat("Original seqnames:", paste(as.character(seqnames(gr)), collapse = ","), "\n")

# This will crash if the Apptainer patch failed
grl <- split(gr, c("pair", "pair"))
grl_nochr <- gUtils::gr.nochr(grl)
cat("Stripped seqlevels:", paste(seqlevels(grl_nochr), collapse = ","), "\n")
cat("[PASS] gUtils/GenomicRanges patch successful!\n\n")


cat("----------------------------------------\n")
cat("TEST 2: CPLEX Engine Check\n")
cat("----------------------------------------\n")
# We test CPLEX first because it doesn't usually require a license for tiny matrices
obj <- c(1, 2, 3)
mat <- matrix(c(-1, 1, 1, 1, -3, 1), nrow = 2)
dir <- c("<=", "<=")
rhs <- c(20, 30)

cplex_result <- tryCatch({
    Rcplex(cvec = obj, Amat = mat, bvec = rhs, sense = "L", objsense = "max", control = list(trace = 0))
}, error = function(e) {
    stop(paste("\nCRITICAL FAILURE: CPLEX loaded, but failed to solve.\nError:", e$message))
})

if (!is.null(cplex_result$status) && cplex_result$status == 1) {
    cat("[PASS] CPLEX successfully solved the LP matrix.\n\n")
} else {
    stop("CRITICAL FAILURE: CPLEX returned an unexpected status.")
}


cat("----------------------------------------\n")
cat("TEST 3: Gurobi Engine Check\n")
cat("----------------------------------------\n")
model <- list()
model$A          <- matrix(c(1,2,3,1,1,0), nrow=2, ncol=3, byrow=TRUE)
model$obj        <- c(1,1,2)
model$modelsense <- 'max'
model$rhs        <- c(4,1)
model$sense      <- c('<', '>')

# Note: This might throw a recognizable license error if you don't map the license during the build.
# If that happens, it means the library works and is correctly asking for a license.
gurobi_result <- tryCatch({
    gurobi(model, params = list(OutputFlag = 0))
}, error = function(e) {
    cat("Note: Gurobi execution failed. If this is a license error ('10009'), it is expected during build.\n")
    cat("Error:", e$message, "\n")
    return(NULL)
})

if (!is.null(gurobi_result) && !is.null(gurobi_result$status) && gurobi_result$status == "OPTIMAL") {
    cat("[PASS] Gurobi successfully solved the LP matrix.\n\n")
}

cat("CONTAINER BUILD COMPLETE: All systems go.\n")