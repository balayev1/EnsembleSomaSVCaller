FROM quay.io/biocontainers/r-sequenza:3.0.0--r42h3342da4_5

LABEL org.opencontainers.image.title="r-sequenza-hg38" \
      org.opencontainers.image.description="Sequenza container with ShixiangWang/copynumber fork for hg38 support" \
      org.opencontainers.image.source="https://github.com/ShixiangWang/copynumber"

ARG COPYNUMBER_REF=2e31d59f558d5e3e44215da2f612707791e3f030

ENV VROOM_CONNECTION_SIZE=33554432

RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("remotes")'

RUN Rscript -e 'remotes::install_github("ShixiangWang/copynumber", ref = Sys.getenv("COPYNUMBER_REF"), upgrade = "never", dependencies = TRUE, build_vignettes = FALSE)'

RUN Rscript -e 'library(copynumber); cat("copynumber:", as.character(packageVersion("copynumber")), "\\n")'
RUN Rscript -e 'library(sequenza); cat("sequenza:", as.character(packageVersion("sequenza")), "\\n")'

COPY test_hg38_extract.R /usr/local/bin/test_hg38_extract.R
RUN chmod 755 /usr/local/bin/test_hg38_extract.R
