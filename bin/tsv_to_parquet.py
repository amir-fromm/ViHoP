#!/usr/bin/env python
# coding: utf-8
import pandas as pd
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser(description="Parallel BLAST summarizer for large files")
    parser.add_argument("--tsv_path", required=True, help="Tsv file name")
    return parser.parse_args()

def convert_blast_tsv_to_parquet(tsv_path):
    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()
    dtypes = {
        "qseqid": "string",
        "sseqid": "string",
        "pident": "Float64",
        "length": "Int64",
        "mismatch": "Int64",
        "gapopen": "Int64",
        "qstart": "Int64",
        "qend": "Int64",
        "sstart": "Int64",
        "send": "Int64",
        "evalue": "Float64",
        "bitscore": "Float64",
    }

    print("Reading TSV to convert...")
    df = pd.read_table(tsv_path, sep="\t", names=cols, dtype=dtypes, on_bad_lines='warn')

    parquet_path = tsv_path.replace('tsv', 'parquet')
    df.to_parquet(parquet_path,engine="fastparquet", compression="snappy", index=False)

def main():
    args = get_args()
    tsv_path = args.tsv_path
    convert_blast_tsv_to_parquet(tsv_path)

if __name__ == "__main__":
    main()
