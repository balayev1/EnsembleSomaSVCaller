#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_FILES = {
    "cnv_tab": "{sample_id}.cnv.tab.gz",
    "snp_tab": "{sample_id}.snp.tab.gz",
    "pscbs_data": "{sample_id}_pscbs_data.txt.gz",
    "sex_file": "{sample_id}_sex.txt",
    "clustered_and_pruned_and_normal": "{sample_id}_clustered_and_pruned_and_normal.txt",
    "ploidy_purity_2d": "{sample_id}_ploidy_purity_2D.txt",
    "breakpoints2": "{sample_id}_breakpoints2.txt",
    "known_segments": "{sample_id}.knownSegments.txt",
}


def resolve_sample_dir(root: Path, sample_id: str) -> Path:
    direct = root / sample_id
    if direct.is_dir():
        return direct

    if root.is_dir() and any((root / pattern.format(sample_id=sample_id)).exists() for pattern in REQUIRED_FILES.values()):
        return root

    for candidate in root.iterdir():
        if candidate.is_dir() and candidate.name == sample_id:
            return candidate

    raise FileNotFoundError(
        f"Could not locate ACEseq sample directory for '{sample_id}' under '{root}'."
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a one-row aceseq_manifest.tsv from a sample ACEseq output directory."
    )
    parser.add_argument("--sample-id", required=True, help="Sample identifier used in the samplesheet.")
    parser.add_argument("--aceseq-root", required=True, help="ACEseq output root for the current run.")
    parser.add_argument("--output", required=True, help="Path to write aceseq_manifest.tsv.")
    args = parser.parse_args()

    sample_id = args.sample_id
    aceseq_root = Path(args.aceseq_root).resolve()
    output = Path(args.output).resolve()

    sample_dir = resolve_sample_dir(aceseq_root, sample_id)

    row = {
        "sample_id": sample_id,
        "aceseq_dir": str(sample_dir),
    }

    missing = []
    for column, pattern in REQUIRED_FILES.items():
        path = sample_dir / pattern.format(sample_id=sample_id)
        if not path.exists():
            missing.append(str(path))
        row[column] = str(path)

    if missing:
        sys.stderr.write("Missing required ACEseq outputs:\n")
        for path in missing:
            sys.stderr.write(f"  - {path}\n")
        return 1

    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sample_id", "aceseq_dir", *REQUIRED_FILES.keys()]

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    print(f"Wrote ACEseq manifest for {sample_id} to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
