#!/usr/bin/env python3
"""Split a master candidate cell list into per-sample TSV files."""

import argparse
from collections import defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import load_candidate_cells


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidates", help="Master candidate cell TSV")
    parser.add_argument("outputs", nargs="+", help="Per-sample output TSV paths")
    args = parser.parse_args()

    grouped: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for barcode, sample in load_candidate_cells(args.candidates):
        grouped[sample].append((barcode, sample))

    for output_path in args.outputs:
        path = Path(output_path)
        sample = path.name.removeprefix("candidate_cells_").removesuffix(".tsv")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as handle:
            for barcode, cell_sample in grouped.get(sample, []):
                handle.write(f"{barcode}\t{cell_sample}\n")


if __name__ == "__main__":
    main()
