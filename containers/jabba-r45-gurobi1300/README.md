# JaBbA + Gurobi 13.0.1 Apptainer image

This directory contains a reproducible Apptainer recipe for a JaBbA image built on:

- `R 4.5`
- `Bioconductor 3.21`
- `JaBbA`
- `slam`
- `Gurobi 13.0.1`

The goal is to replace the older `mskilab/jabba` runtime with an image that can
load the `gurobi` R package directly and carry the JaBbA compatibility fixes in
the container startup, so the Nextflow module does not need a runtime monkey
patch.

## Why this image stays on R 4.5 / Bioconductor 3.21

Gurobi 13.0.1 ships a Linux R package for the R 4.5 series, and the local
installer bundled for this build contains:

```text
/opt/gurobi1301/linux64/R/gurobi_13.0-1_R_4.5.0.tar.gz
```

So this image keeps the matching `R 4.5` / `Bioconductor 3.21` base and applies
the JaBbA/gGnome compatibility shim through `/usr/local/lib/R/etc/Rprofile.site`
inside the container.

## Important note

This repo intentionally does **not** include:

- any `gurobi.lic` file
- any Gurobi installer tarball

You must provide your own local Gurobi installer at build time and pass the
license at runtime.

If you publish the resulting image to a public registry such as Quay, make sure
your Gurobi license and redistribution terms allow that.

## Required local file before build

Place the Gurobi Linux installer tarball in this directory as:

```text
containers/jabba-r45-gurobi1300/gurobi13.0.1_linux64.tar.gz
```

The recipe expects that tarball to unpack into:

```text
/opt/gurobi1300/linux64
```

and contain:

```text
/opt/gurobi1300/linux64/R/gurobi_13.0-1_R_4.5.0.tar.gz
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

The image includes a simple R smoke test:

```bash
apptainer exec \
  --env GRB_LICENSE_FILE=/path/to/gurobi.lic \
  jabba-r45-gurobi1300_0.2.0.sif \
  Rscript /usr/local/bin/test_jabba_gurobi.R
```

It checks that `JaBbA`, `slam`, and `gurobi` all load successfully, and that the
container-scoped compatibility shim restores the `GRanges`/`GRangesList` method
dispatch JaBbA needs.

## Suggested pipeline follow-up

After publishing the image, the pipeline can be updated to use:

```text
oras://quay.io/<quay-user>/jabba-r45-gurobi1300:0.2.0
```

for the JaBbA process instead of the current hard-coded image.
