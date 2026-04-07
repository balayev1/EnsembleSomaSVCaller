#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = [
    "sample_id",
    "aceseq_dir",
    "cnv_tab",
    "snp_tab",
    "pscbs_data",
    "sex_file",
    "clustered_and_pruned_and_normal",
    "ploidy_purity_2d",
    "breakpoints2",
    "known_segments",
]


def load_samplesheet(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "sample" not in reader.fieldnames:
            raise ValueError(f"Samplesheet is missing required 'sample' column: {path}")
        rows = list(reader)
    return [row["sample"] for row in rows]


def load_manifest(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing_columns = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
        if missing_columns:
            raise ValueError(
                f"ACEseq manifest is missing required columns: {', '.join(missing_columns)}"
            )

        manifest = {}
        for row in reader:
            sample_id = row["sample_id"]
            if sample_id in manifest:
                raise ValueError(f"Duplicate sample_id in ACEseq manifest: {sample_id}")
            manifest[sample_id] = row

    return manifest


def validate_paths(sample_id: str, row: dict):
    missing = []
    for column in REQUIRED_COLUMNS[1:]:
        value = row.get(column, "").strip()
        if not value:
            missing.append(f"{column} is empty")
            continue
        if not Path(value).exists():
            missing.append(f"{column} missing: {value}")
    if missing:
        raise FileNotFoundError(
            f"Sample {sample_id} has invalid ACEseq artifacts:\n  - " + "\n  - ".join(missing)
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate that the ACEseq manifest contains all required files for the samplesheet."
    )
    parser.add_argument("--samplesheet", required=True, help="SomaSV CSV samplesheet.")
    parser.add_argument("--manifest", required=True, help="ACEseq manifest TSV.")
    args = parser.parse_args()

    samplesheet = Path(args.samplesheet).resolve()
    manifest_path = Path(args.manifest).resolve()

    samples = load_samplesheet(samplesheet)
    manifest = load_manifest(manifest_path)

    missing_samples = [sample for sample in samples if sample not in manifest]
    if missing_samples:
        sys.stderr.write(
            "ACEseq manifest is missing samples required by the SomaSV samplesheet:\n"
        )
        for sample in missing_samples:
            sys.stderr.write(f"  - {sample}\n")
        return 1

    for sample in samples:
        validate_paths(sample, manifest[sample])

    print(
        f"Validated ACEseq manifest {manifest_path} for {len(samples)} sample(s) from {samplesheet}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
