#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from argparse import ArgumentParser
from multiprocessing import Pool, cpu_count


def get_args():
    parser = ArgumentParser(description="Parallel BLAST summarizer for large files")
    parser.add_argument("--path", required=True, help="Path to all input and output files")
    parser.add_argument("--metapr2_file", required=True, help="MetaPR2 BLAST file name")
    parser.add_argument("--pr2_file", required=True, help="PR2 BLAST file name")
    parser.add_argument("--replacement_file", required=True, help="Replacement table (TSV)")
    parser.add_argument("--identity", type=float, default=99, help="Minimum percent identity")
    parser.add_argument("--length", type=float, default=100, help="Minimum alignment length")
    parser.add_argument("--level", default="Annotation", help="Taxonomic level to summarize")
    parser.add_argument("--cpus_available", type=int, default=1, help="Number of cpus available")
    return parser.parse_args()


def load_data_chunked(file_path, db, chunk_size=5_000_000):
    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(" ")
    chunks = pd.read_table(file_path, sep="\t", names=cols, chunksize=chunk_size)
    for chunk in chunks:
        chunk['cell'] = chunk['qseqid'].str.split("i\\d+_", expand=True)[1]
        if db == "PR2":
            chunk = chunk[chunk['sseqid'].str.contains("18S_rRNA")]
            chunk[["sseqid", "rrna_type", "rrna_loc", "clone", "Domain", "Clade", "Supergroup",
                   "Division", "Class", "Order", "Family", "Genus", "Species"]] = chunk['sseqid'].str.split("|", expand=True)
        elif db == "metaPR2":
            chunk[["sseqid", "Clade", "Supergroup", "Division", "Class", "Order",
                   "Family", "Genus", "Species"]] = chunk['sseqid'].str.split("|", expand=True)
        yield chunk


def replace_worker(args):
    df_top, level, old_name, new_name = args
    group = set(df_top[df_top[level] == old_name]['Family'])
    return {x: new_name for x in group}


def compute_super_dict(all_data, replacement_df):
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


def process_chunk(chunk, super_dict, identity, length, level):
    chunk['Annotation'] = chunk['Family'].map(super_dict)
    chunk['barcode'] = chunk['cell'].str.split(".", expand=True)[0]
    chunk = subset_results(chunk, identity, length)
    chunk = chunk.drop_duplicates(subset=["qseqid", level]).drop_duplicates('qseqid', keep=False)
    chunk = get_chimeric_cells(chunk, level)
    return chunk


def parallel_process_all_chunks(file_path, db, super_dict, identity, length, level,cpus):
    chunk_gen = load_data_chunked(file_path, db)
    pool = Pool(processes=min(cpu_count(), cpus))
    results = []
    for chunk in chunk_gen:
        results.append(pool.apply_async(process_chunk, (chunk, super_dict, identity, length, level)))
    pool.close()
    pool.join()
    return pd.concat([r.get() for r in results])


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
    cpus = args.cpus_available
    replacement_df = pd.read_table(path + args.replacement_file, sep="\t")
    mini_meta = next(load_data_chunked(path + args.metapr2_file, "metaPR2", 500_000))
    mini_pr2 = next(load_data_chunked(path + args.pr2_file, "PR2", 500_000))
    all_data = pd.concat([mini_meta, mini_pr2])
    super_dict = compute_super_dict(all_data, replacement_df)

    print("Processing MetaPR2...")
    metapr2_data = parallel_process_all_chunks(path + args.metapr2_file, "metaPR2", super_dict, identity, length, level, cpus)
    print("Processing PR2...")
    pr2_data = parallel_process_all_chunks(path + args.pr2_file, "PR2", super_dict, identity, length, level,cpus)

    print("Annotating cells...")
    result_df = annotation_data(pr2_data, metapr2_data, level)
    result_df.to_csv(path + f"blast_summary.{level}.csv", index=False)


if __name__ == "__main__":
    main()
