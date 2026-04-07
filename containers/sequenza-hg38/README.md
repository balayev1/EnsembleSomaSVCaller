# Sequenza hg38 image

This directory contains the reproducible recipe for the custom Sequenza image
used by `SEQUENZA_RUN`.

It is based on `r-sequenza:3.0.0` and installs the
`ShixiangWang/copynumber` fork pinned to:

- `2e31d59f558d5e3e44215da2f612707791e3f030`

The built `.sif` is intentionally not committed. Rebuild it from the files in
this directory on any new system.

The currently published image is:

- `oras://quay.io/balay011/r-sequenza-hg38:3.0.0-copynumber-2e31d59`

## Local rebuild option

If you want a local `.sif` copy on another system, build it at:

```text
containers/sequenza-hg38/r-sequenza-hg38_3.0.0-copynumber-2e31d59.sif
```

After cloning the repo elsewhere, you can build or copy the image into this same
directory with that filename.

## Build with Apptainer on HPC

If the HPC allows unprivileged builds with `--fakeroot`, run:

```bash
cd containers/sequenza-hg38
apptainer build --fakeroot \
  r-sequenza-hg38_3.0.0-copynumber-2e31d59.sif \
  apptainer.def
```

## Build with Docker elsewhere, then transfer

If the target HPC does not allow `apptainer build --fakeroot`, build the image
on a machine with Docker:

```bash
cd containers/sequenza-hg38
docker build -t r-sequenza-hg38:3.0.0-copynumber-2e31d59 .
```

Then either:

1. Convert it to a `.sif` on a system with Apptainer/Singularity.
2. Push it to a registry and pull it back as a `.sif` on the HPC.

Example pull from a registry-backed image:

```bash
cd containers/sequenza-hg38
singularity pull r-sequenza-hg38_3.0.0-copynumber-2e31d59.sif \
  docker://ghcr.io/<user>/r-sequenza-hg38:3.0.0-copynumber-2e31d59
```

## Smoke test

After building the image, run:

```bash
cd containers/sequenza-hg38
singularity exec r-sequenza-hg38_3.0.0-copynumber-2e31d59.sif \
  Rscript /usr/local/bin/test_hg38_extract.R /path/to/CB822.merged.seqz.gz
```

## Notes

- The pipeline currently uses the published Quay ORAS image in
  `SEQUENZA_RUN`.
- The files in this directory are kept so the image can be rebuilt or updated
  reproducibly later.
- The pipeline itself sets `VROOM_CONNECTION_SIZE` inside the R task script, so
  the runtime behavior does not depend only on the container environment.
