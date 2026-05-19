#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

MIN_CALLERS_FOR_CONSENSUS = 2


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge purity/ploidy calls from ASCAT, FACETS, Sequenza, and ACEseq."
    )
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--ascat", required=True, help="ASCAT *.purityploidy.txt file")
    parser.add_argument("--facets", required=True, help="FACETS *.vcf.gz file")
    parser.add_argument(
        "--sequenza",
        required=True,
        help="Sequenza *_alternative_solutions.txt file",
    )
    parser.add_argument(
        "--aceseq",
        required=False,
        help="Optional ACEseq *_ploidy_purity_2D.txt file",
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=0.1,
        help="DBSCAN epsilon distance threshold",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=2,
        help="DBSCAN minimum number of points to define a dense region",
    )
    parser.add_argument(
        "--candidates-out",
        required=True,
        help="Output TSV with candidate calls and cluster labels",
    )
    parser.add_argument(
        "--consensus-out",
        required=True,
        help="Output TSV with one consensus row for the sample",
    )
    return parser.parse_args()


def require_float(value, description):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Could not parse numeric value for {description}: {value!r}") from exc

    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite numeric value for {description}: {value!r}")

    return parsed


def is_missing_value(value):
    if value is None:
        return True

    text = str(value).strip()
    return text == "" or text.upper() in {"NA", "NULL", "NAN", "."}


def read_tsv_rows(path):
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            cleaned = {key: (value.strip() if value is not None else "") for key, value in row.items()}
            if any(cleaned.values()):
                rows.append(cleaned)
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return rows


def make_candidate(sample_id, caller, base_caller, option_rank, purity, ploidy, source_file):
    if is_missing_value(purity) or is_missing_value(ploidy):
        return None

    return {
        "sample_id": sample_id,
        "caller": caller,
        "base_caller": base_caller,
        "option_rank": option_rank,
        "purity": require_float(purity, f"{caller} purity"),
        "ploidy": require_float(ploidy, f"{caller} ploidy"),
        "source_file": str(Path(source_file).resolve()),
    }


def parse_ascat(sample_id, path):
    row = read_tsv_rows(path)[0]
    candidate = make_candidate(
        sample_id=sample_id,
        caller="ASCAT",
        base_caller="ASCAT",
        option_rank=1,
        purity=row["AberrantCellFraction"],
        ploidy=row["Ploidy"],
        source_file=path,
    )
    return [candidate] if candidate is not None else []


def parse_facets(sample_id, path):
    purity = None
    ploidy = None

    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("##purity="):
                purity = line.split("=", 1)[1].strip()
            elif line.startswith("##ploidy="):
                ploidy = line.split("=", 1)[1].strip()
            elif line.startswith("#CHROM"):
                break

    if purity is None or ploidy is None:
        raise ValueError(f"Could not find ##purity or ##ploidy headers in {path}")

    candidate = make_candidate(
        sample_id=sample_id,
        caller="FACETS",
        base_caller="FACETS",
        option_rank=1,
        purity=purity,
        ploidy=ploidy,
        source_file=path,
    )
    return [candidate] if candidate is not None else []


def parse_sequenza(sample_id, path):
    row = read_tsv_rows(path)[0]
    candidate = make_candidate(
        sample_id=sample_id,
        caller="Sequenza",
        base_caller="Sequenza",
        option_rank=1,
        purity=row["cellularity"],
        ploidy=row["ploidy"],
        source_file=path,
    )
    return [candidate] if candidate is not None else []


def parse_aceseq(sample_id, path):
    rows = read_tsv_rows(path)
    candidates = []
    for idx, row in enumerate(rows, start=1):
        candidate = make_candidate(
            sample_id=sample_id,
            caller=f"ACEseq_Option{idx}",
            base_caller="ACEseq",
            option_rank=idx,
            purity=row["tcc"],
            ploidy=row["ploidy_factor"],
            source_file=path,
        )
        if candidate is not None:
            candidates.append(candidate)
    return candidates


def euclidean_distance(point_a, point_b):
    return math.sqrt((point_a[0] - point_b[0]) ** 2 + (point_a[1] - point_b[1]) ** 2)


def region_query(points, point_index, eps):
    center = points[point_index]
    return [idx for idx, point in enumerate(points) if euclidean_distance(center, point) <= eps]


def dbscan(points, eps, min_samples):
    labels = [None] * len(points)
    cluster_id = 0

    for point_index in range(len(points)):
        if labels[point_index] is not None:
            continue

        neighbors = region_query(points, point_index, eps)
        if len(neighbors) < min_samples:
            labels[point_index] = -1
            continue

        labels[point_index] = cluster_id
        seeds = [idx for idx in neighbors if idx != point_index]
        seed_cursor = 0

        while seed_cursor < len(seeds):
            neighbor_index = seeds[seed_cursor]

            if labels[neighbor_index] == -1:
                labels[neighbor_index] = cluster_id

            if labels[neighbor_index] is not None:
                seed_cursor += 1
                continue

            labels[neighbor_index] = cluster_id
            neighbor_neighbors = region_query(points, neighbor_index, eps)
            if len(neighbor_neighbors) >= min_samples:
                for idx in neighbor_neighbors:
                    if idx not in seeds:
                        seeds.append(idx)
            seed_cursor += 1

        cluster_id += 1

    return labels


def mean_distance_to_centroid(records):
    centroid = (
        statistics.mean(record["purity"] for record in records),
        statistics.mean(record["ploidy"] for record in records),
    )
    return statistics.mean(
        euclidean_distance((record["purity"], record["ploidy"]), centroid) for record in records
    )


def unique_callers_in_order(records):
    seen = set()
    callers = []
    for record in records:
        if record["base_caller"] not in seen:
            seen.add(record["base_caller"])
            callers.append(record["base_caller"])
    return callers


def count_missing_base_callers(candidates, expected_callers):
    usable_callers = set(unique_callers_in_order(candidates))
    return sum(1 for caller in expected_callers if caller not in usable_callers)


def assign_clusters(candidates, eps, min_samples):
    points = [(candidate["purity"], candidate["ploidy"]) for candidate in candidates]
    labels = dbscan(points, eps=eps, min_samples=min_samples)

    for candidate, label in zip(candidates, labels):
        candidate["cluster_label"] = label
        candidate["in_consensus"] = False

    clusters = defaultdict(list)
    for candidate in candidates:
        if candidate["cluster_label"] >= 0:
            clusters[candidate["cluster_label"]].append(candidate)

    valid_clusters = {
        label: records
        for label, records in clusters.items()
        if len(set(record["base_caller"] for record in records)) >= 2
    }

    if not valid_clusters:
        return None

    def cluster_rank(item):
        label, records = item
        return (
            len(set(record["base_caller"] for record in records)),
            len(records),
            -mean_distance_to_centroid(records),
            -label,
        )

    best_label, best_records = max(valid_clusters.items(), key=cluster_rank)
    for candidate in candidates:
        candidate["in_consensus"] = candidate["cluster_label"] == best_label

    return best_label, best_records


def write_candidates(path, candidates):
    fieldnames = [
        "sample_id",
        "caller",
        "base_caller",
        "option_rank",
        "purity",
        "ploidy",
        "cluster_label",
        "in_consensus",
        "source_file",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for candidate in candidates:
            writer.writerow(candidate)


def write_consensus(path, sample_id, consensus, eps, min_samples):
    fieldnames = [
        "sample_id",
        "status",
        "consensus_cluster_label",
        "purity",
        "ploidy",
        "supporting_callers",
        "supporting_points",
        "n_supporting_callers",
        "n_supporting_points",
        "eps",
        "min_samples",
    ]

    row = {
        "sample_id": sample_id,
        "status": "no_consensus",
        "consensus_cluster_label": "NA",
        "purity": "NA",
        "ploidy": "NA",
        "supporting_callers": "",
        "supporting_points": "",
        "n_supporting_callers": 0,
        "n_supporting_points": 0,
        "eps": eps,
        "min_samples": min_samples,
    }

    if consensus is not None:
        label, records = consensus
        supporting_callers = unique_callers_in_order(records)
        row.update(
            {
                "status": "consensus",
                "consensus_cluster_label": label,
                "purity": statistics.median(record["purity"] for record in records),
                "ploidy": statistics.median(record["ploidy"] for record in records),
                "supporting_callers": ",".join(supporting_callers),
                "supporting_points": ",".join(record["caller"] for record in records),
                "n_supporting_callers": len(supporting_callers),
                "n_supporting_points": len(records),
            }
        )

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def main():
    args = parse_args()
    expected_callers = ["ASCAT", "FACETS", "Sequenza"]
    if args.aceseq:
        expected_callers.append("ACEseq")

    candidates = []
    candidates.extend(parse_ascat(args.sample_id, args.ascat))
    candidates.extend(parse_facets(args.sample_id, args.facets))
    candidates.extend(parse_sequenza(args.sample_id, args.sequenza))
    if args.aceseq:
        candidates.extend(parse_aceseq(args.sample_id, args.aceseq))

    missing_base_callers = count_missing_base_callers(candidates, expected_callers)
    consensus = assign_clusters(candidates, eps=args.eps, min_samples=args.min_samples)
    if missing_base_callers >= len(expected_callers) - MIN_CALLERS_FOR_CONSENSUS + 1:
        consensus = None

    write_candidates(args.candidates_out, candidates)
    write_consensus(args.consensus_out, args.sample_id, consensus, args.eps, args.min_samples)

    print(f"Collected {len(candidates)} purity/ploidy candidates for sample {args.sample_id}.")
    if not candidates:
        print("No usable purity/ploidy candidates remained after filtering missing values.")
    elif missing_base_callers >= len(expected_callers) - MIN_CALLERS_FOR_CONSENSUS + 1:
        threshold = len(expected_callers) - MIN_CALLERS_FOR_CONSENSUS + 1
        print(
            f"No consensus cluster was reported because {threshold} or more expected callers yielded missing purity/ploidy values."
        )
    elif consensus is None:
        print("No consensus cluster was found across at least two distinct callers.")
    else:
        _, records = consensus
        purity = statistics.median(record["purity"] for record in records)
        ploidy = statistics.median(record["ploidy"] for record in records)
        callers = ", ".join(unique_callers_in_order(records))
        print(f"Consensus purity={purity:.6f} ploidy={ploidy:.6f} supported by {callers}.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
