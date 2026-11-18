#!/usr/bin/env python3
"""Assign an NCBI division to the best DIAMOND alignment of each query.

The alignment table is DIAMOND tabular output, 12 columns, as written by:

    diamond blastp --outfmt 6 qseqid sseqid pident length mismatch gapopen \\
        qstart qend sstart send evalue bitscore

The same 12 columns are what BLAST and lastal -f BlastTab write. Extra
columns after the bit score are ignored. A header line (starting with #,
or with the field name qseqid) is skipped. For each query the rows with
the highest bit score are kept. Those subject accessions are looked up in
the RefSeq taxid catalog with one sequential read, then mapped to Bacteria,
Eukaryotes, Phages, Viruses, or Others.

A subject missing from the catalog is labeled Viruses. A query whose tied
best hits fall in more than one division is removed.
"""

import argparse
import sys
from collections import OrderedDict

COLUMNS = [
    "query id",
    "subject id",
    "% identity",
    "alignment length",
    "mismatches",
    "gap opens",
    "q. start",
    "q. end",
    "s. start",
    "s. end",
    "evalue",
    "bit score",
    "rank",
    "Domain",
]

DIVISIONS = {
    0: "Bacteria",
    1: "Eukaryotes",
    2: "Eukaryotes",
    3: "Phages",
    4: "Eukaryotes",
    5: "Eukaryotes",
    6: "Eukaryotes",
    7: "Others",
    8: "Others",
    9: "Viruses",
    10: "Eukaryotes",
    11: "Others",
}


def best_hits(path, excluded_queries):
    """Keep every alignment that ties the highest bit score of its query."""
    best = OrderedDict()
    alignments = 0
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[0] in ("qseqid", "query", "query id"):
                continue
            query = fields[0]
            if query in excluded_queries:
                continue
            score = fields[11]
            try:
                score_value = float(score)
            except ValueError:
                continue
            alignments += 1
            current = best.get(query)
            row = fields[:12]
            if current is None or score_value > current[0]:
                best[query] = (score_value, [row])
            elif score_value == current[0]:
                current[1].append(row)
    print(f"Alignments read: {alignments:,}; queries with a hit: {len(best):,}", flush=True)
    return best


def load_nodes(path):
    nodes = {}
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                nodes[parts[0]] = int(parts[1])
            except ValueError:
                continue
    return nodes


def load_merged(path):
    merged = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                merged[parts[0]] = parts[1]
    return merged


def division_name(taxid, nodes, merged):
    """Return the division, or None when the taxid cannot be labeled."""
    def from_nodes(identifier):
        division = nodes.get(identifier)
        if division is None:
            return None
        return DIVISIONS.get(division)

    found = from_nodes(taxid)
    if found is not None:
        return found
    new_taxid = merged.get(taxid)
    if new_taxid is not None:
        print(f"Merged taxid {taxid} : {new_taxid}", flush=True)
        found = from_nodes(new_taxid)
        if found is not None:
            return found
        print(f"{new_taxid} taxid not in nodes dictionary", flush=True)
        return None
    print(f"Unknown taxid {taxid}", flush=True)
    return "Unknown taxid"


def lookup_accessions(catalog_path, needed, nodes, merged):
    """Map subject accessions to divisions by scanning the catalog once."""
    found = {}
    seen = 0
    with open(catalog_path) as handle:
        for line in handle:
            seen += 1
            if seen % 50000000 == 0:
                print(f"Catalog lines scanned: {seen:,}; accessions matched: {len(found):,}", flush=True)
            parts = line.split()
            if len(parts) < 2:
                continue
            accession = parts[1]
            if accession not in needed:
                continue
            taxid = parts[0]
            label = division_name(taxid, nodes, merged)
            if label is not None:
                found[accession] = label
            if len(found) == len(needed):
                break
    print(f"Catalog lines scanned: {seen:,}; accessions matched: {len(found):,} of {len(needed):,}", flush=True)
    return found


def rows_to_write(best, accession_division):
    """One row per query, unless tied best hits disagree on the division."""
    written = []
    dropped = 0
    for _query, (_score, rows) in best.items():
        labeled = []
        for row in rows:
            domain = accession_division.get(row[1], "Viruses")
            labeled.append((row, domain))
        domains = {domain for _row, domain in labeled}
        if len(domains) > 1:
            dropped += 1
            continue
        row, domain = labeled[0]
        written.append(row + ["1", domain])
    print(f"Queries removed because tied hits have different divisions: {dropped:,}", flush=True)
    return written


def write_table(path, rows):
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def read_id_list(path):
    with open(path) as handle:
        return {line.strip() for line in handle if line.strip()}


def main():
    parser = argparse.ArgumentParser(
        description="Label the best DIAMOND alignment of each query with an NCBI division."
    )
    parser.add_argument(
        "--alignment",
        "--diamond_path",
        "--last_path",
        dest="alignment",
        required=True,
        help="DIAMOND outfmt 6 table (12 columns: qseqid through bitscore). "
        "--diamond_path and --last_path are the same argument.",
    )
    parser.add_argument("--output_path", default="summary_diamond.tsv")
    parser.add_argument("--chunks", type=float, default=100000, help="Accepted for compatibility; the file is streamed")
    parser.add_argument("--size_table", type=float, default=0, help="Accepted for compatibility; the file is streamed")
    parser.add_argument("--false_positives", default="-1", help="Query IDs to remove before labeling")
    parser.add_argument("--catalog", required=True, help="RefSeq taxid catalog: taxid, space, accession")
    parser.add_argument("--nodes", required=True, help="taxid and NCBI division, whitespace-separated")
    parser.add_argument("--merged", required=True, help="old taxid, tab, new taxid")
    args = parser.parse_args()

    excluded = set()
    if args.false_positives not in ("-1", "", None):
        excluded = read_id_list(args.false_positives)
        print(f"Excluding {len(excluded):,} query IDs", flush=True)

    best = best_hits(args.alignment, excluded)
    if not best:
        write_table(args.output_path, [])
        print(f"Wrote 0 rows -> {args.output_path}", flush=True)
        return

    needed = {row[1] for _score, rows in best.values() for row in rows}
    print(f"Unique subjects to look up: {len(needed):,}", flush=True)
    nodes = load_nodes(args.nodes)
    merged = load_merged(args.merged)
    accession_division = lookup_accessions(args.catalog, needed, nodes, merged)
    rows = rows_to_write(best, accession_division)
    write_table(args.output_path, rows)
    print(f"Wrote {len(rows):,} rows -> {args.output_path}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
