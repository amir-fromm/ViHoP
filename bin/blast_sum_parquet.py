#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from argparse import ArgumentParser
from multiprocessing import Pool, cpu_count
from fastparquet import ParquetFile
import sys

def get_args():
    parser = ArgumentParser(description="Parallel BLAST summarizer for large files")
    parser.add_argument("--path", required=True, help="Path to all input and output files")
    parser.add_argument("--metapr2_file", required=True, help="MetaPR2 Parquet file name")
    parser.add_argument("--pr2_file", required=True, help="PR2 Parquet file name")
    parser.add_argument("--replacement_file", required=True, help="Replacement table (TSV)")
    parser.add_argument("--identity", type=float, default=99, help="Minimum percent identity")
    parser.add_argument("--length", type=float, default=100, help="Minimum alignment length")
    parser.add_argument("--level", default="Annotation", help="Taxonomic level to summarize")
    return parser.parse_args()


def read_parquet_file(parquet_file):
    pf = ParquetFile(parquet_file)
    df = pf.to_pandas()
    df['sseqid'] = df['sseqid'].astype(str)
    df['cell'] = df['qseqid'].str.extract(r"i\d+_(.+):.+$")[0]
    return df

def split_sseqid(df,db):
    if db == "PR2":
        df = df[df['sseqid'].str.contains("18S_rRNA", na=False)]
        split_cols = df['sseqid'].str.split("|", expand=True)
        if split_cols.shape[1] >= 13:
            df[["sseqid", "rrna_type", "rrna_loc", "clone", "Domain", "Clade", "Supergroup",
                "Division", "Class", "Order", "Family", "Genus", "Species"]] = split_cols.iloc[:, :14]
    elif db == "metaPR2":
        split_cols = df['sseqid'].str.split("|", expand=True)
        if split_cols.shape[1] >= 9:
            df[["sseqid", "Clade", "Supergroup", "Division", "Class", "Order",
                "Family", "Genus", "Species"]] = split_cols.iloc[:, :9]
    return df


def replace_worker(args):
    df_top, level, old_name, new_name = args
    if level not in df_top.columns:
        return {}
    group = set(df_top[df_top[level] == old_name]['Family'])
    return {x: new_name for x in group}


def compute_super_dict(all_data, replacement_df):
    print (all_data.columns)
    print("Building replacement dictionary...")
    tasks = [(all_data, row['level'], row['old_name'], row['new_name']) for _, row in replacement_df.iterrows()]
    results = [replace_worker(task) for task in tasks]
    super_dict = {k: v for d in results for k, v in d.items()}
    rest = {x: "Other eukaryotes" for x in all_data.Family.unique() if x not in super_dict}
    super_dict.update(rest)
    return super_dict


def subset_results(df, identity, length):
    df = df[df['length'] >= length]
    df = df.sort_values(by=['qseqid', 'bitscore'], ascending=False)
    df['rank'] = df.groupby('qseqid')['bitscore'].rank("min", ascending=False)
    df = df[df['rank'] == 1]
    df = df[df['pident'] >= identity]
    return df


def find_chimeras(df, level):
    df = df.drop_duplicates(subset=["qseqid", level])
    df = df.sort_values(by=['qseqid', 'qstart'])
    result = set()
    curr_id, curr_start, curr_end = None, None, None
    ranges = set()
    for _, row in df.iterrows():
        gid, start, end = row['qseqid'], row['qstart'], row['qend']
        if gid != curr_id:
            if curr_id and len(ranges) >= 2:
                result.add(curr_id)
            curr_id, curr_start, curr_end = gid, start, end
            ranges = {f"{start}-{end}"}
        else:
            ranges.add(f"{start}-{end}")
            curr_end = max(curr_end, end)
    return result


def get_chimeric_cells(df, level):
    chimeras = find_chimeras(df, level)
    return df[~df['cell'].isin(df[df['qseqid'].isin(chimeras)]['cell'])]


def process_chunk(chunk, replacement_df, identity, length, level, db):
    chunk = subset_results(chunk, identity, length)
    chunk = split_sseqid(chunk, db)
    super_dict = compute_super_dict(chunk, replacement_df)
    chunk['Annotation'] = chunk['Family'].map(super_dict)
    chunk['barcode'] = chunk['cell'].str.split(".", expand=True)[0]
    chunk = chunk.drop_duplicates(subset=["qseqid", level]).drop_duplicates('qseqid', keep=False)
    chunk = get_chimeric_cells(chunk, level)
    return chunk


def parallel_process_df(full_df, replacement_df, identity, length, level, db):
    num_chunks = min(cpu_count(), 8)
    chunks = np.array_split(full_df, num_chunks)
    args = [(chunk.copy(), replacement_df, identity, length, level,db) for chunk in chunks]
    with Pool(num_chunks) as pool:
        results = pool.starmap(process_chunk, args)
    return pd.concat(results)


def annotation_data(pr2_data, metapr2_data, level):
    sub1 = pr2_data[["qseqid", "cell", level, "length"]].copy()
    sub2 = metapr2_data[["qseqid", "cell", level, "length"]].copy()
    sub1['cell'] = sub1['cell'].str.strip()
    sub2['cell'] = sub2['cell'].str.strip()
    merged = pd.concat([sub1, sub2])
    merged = merged.sort_values(by='length', ascending=False).drop_duplicates(subset=["qseqid", level])
    merged['percent'] = merged['length'] / merged.groupby('cell')['length'].transform('sum')

    return merged.groupby(['cell', level]).agg({'percent': 'sum', 'length': 'sum'}).reset_index()


def main():
    args = get_args()
    path = args.path.rstrip("/") + "/"
    level = args.level
    identity = args.identity
    length = args.length

    replacement_df = pd.read_table(path + args.replacement_file, sep="\t")

    print("Reading MetaPR2...")
    meta_df = read_parquet_file(path + args.metapr2_file)
    print("Reading PR2...")
    pr2_df = read_parquet_file(path + args.pr2_file)
  #  all_data = pd.concat([meta_df, pr2_df])

    

    print("Processing PR2...")
    pr2_data = parallel_process_df(pr2_df,replacement_df, identity, length, level,'PR2')
    print("Processing MetaPR2...")
    metapr2_data = parallel_process_df(meta_df,replacement_df, identity, length, level,'metaPR2')
    print("Annotating cells...")
    print(pr2_data.shape,metapr2_data.shape)
    result_df = annotation_data(pr2_data, metapr2_data, level)
    print(result_df.shape)
    result_df.to_csv(path + f"blast_summary.{level}.csv", index=False)


if __name__ == "__main__":
    import numpy as np
    main()
