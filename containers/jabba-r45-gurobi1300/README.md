# JaBbA + Gurobi + CPLEX Apptainer image

This directory contains a reproducible Apptainer recipe for a JaBbA image built on:

- `R 4.5`
- `Bioconductor 3.21`
- `JaBbA`
- `slam`
- `Gurobi 13.0.1`
- `IBM CPLEX` with the `Rcplex` R package

The goal is to provide one JaBbA runtime that supports both optimization paths used by this repository:

- `gurobi_jabba=TRUE` via the `gurobi` R package
- `gurobi_jabba=FALSE` via JaBbA's `Rcplex` fallback path

## Why this image stays on R 4.5 / Bioconductor 3.21

Gurobi 13.0.1 ships a Linux R package for the R 4.5 series, so this image keeps
the matching `R 4.5` / `Bioconductor 3.21` base and applies the JaBbA/gGnome
compatibility shim through `/usr/local/lib/R/etc/Rprofile.site` inside the
container.

## Important note

This repo intentionally does **not** include:

- any `gurobi.lic` file
- any Gurobi installer tarball
- any IBM CPLEX installer `.bin`

You must provide your own local solver installers at build time, and your own
Gurobi license at runtime when Gurobi is enabled.

If you publish the resulting image to a public registry such as Quay, make sure
your Gurobi and IBM CPLEX redistribution terms allow that.

## Required local files before build

Place these files in this directory before building:

```text
containers/jabba-r45-gurobi1300/gurobi13.0.1_linux64.tar.gz
containers/jabba-r45-gurobi1300/cos_installer_preview-22.1.2.R4-M0N96ML-linux-x86-64.bin
```

The recipe expects:

- the Gurobi tarball to unpack under `/opt/gurobi1300/linux64`
- the IBM CPLEX `.bin` installer to run in silent mode and install under `/opt/cplex2211`

In particular, the CPLEX layout must provide:

```text
/opt/cplex2211/cplex/include
/opt/cplex2211/cplex/lib/x86-64_linux/static_pic
```

That is what the `Rcplex` installation step uses.

The silent install flow used in the Apptainer recipe is:

```bash
HOME=/tmp/cplex_home \
/opt/cplex_installer.bin \
  -i silent \
  -tempdir /tmp/cplex_installer_tmp \
  -DUSER_INSTALL_DIR=/opt/cplex2211 \
  -DLICENSE_ACCEPTED=TRUE
```

## Build on HPC with Apptainer

If your HPC allows `--fakeroot`, run:

```bash
cd containers/jabba-r45-gurobi1300
apptainer build --fakeroot \
  jabba-r45-gurobi1300_0.2.0.sif \
  apptainer.def
```

If `--fakeroot` is not available, the image will need to be built on another
system with Apptainer root privileges, then copied back to the HPC.

## Push the built SIF to Quay via ORAS

After logging in with Apptainer:

```bash
apptainer registry login quay.io
apptainer push \
  jabba-r45-gurobi1300_0.2.0.sif \
  oras://quay.io/<quay-user>/jabba-r45-gurobi1300:0.2.0
```

## Pull back from Quay

```bash
apptainer pull \
  jabba-r45-gurobi1300_0.2.0.sif \
  oras://quay.io/<quay-user>/jabba-r45-gurobi1300:0.2.0
```

## Smoke test

The image now includes a combined optimizer smoke test:

```bash
apptainer exec \
  --env GRB_LICENSE_FILE=/path/to/gurobi.lic \
  jabba-r45-gurobi1300_0.2.0.sif \
  Rscript /usr/local/bin/test_jabba_optimizers.R
```

It checks that these packages all load successfully:

- `JaBbA`
- `gurobi`
- `Rcplex`
- `slam`

If you want to validate only the non-Gurobi fallback stack, you can still run
the same script without planning to use Gurobi at runtime, but the package load
test itself still requires the Gurobi R package to be present in the image.

## Suggested pipeline follow-up

After publishing the image, the pipeline can continue using:

```text
oras://quay.io/<quay-user>/jabba-r45-gurobi1300:0.2.0
```

for the JaBbA process, but now that image should be able to support both
`gurobi_jabba=TRUE` and `gurobi_jabba=FALSE`, assuming the bundled CPLEX tree is
compatible with `Rcplex`.
