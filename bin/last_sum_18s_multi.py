#!/usr/bin/env python
# coding: utf-8

import glob
import pandas as pd
from tqdm import tqdm
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor

# Load all LAST output tables, combine them, and get only the best hit for each read
def load_last(last_path):
    print("Processing LAST tables")
    flist = glob.glob(last_path + "*.last")

    df_list = []
  #  with ProcessPoolExecutor() as executor:
  #      results = list(tqdm(executor.map(process_file, flist, [last_path] * len(flist)), total=len(flist)))

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for result in tqdm(executor.map(process_file, flist, [last_path] * len(flist)), total=len(flist)):
            if result is not None:
                results.append(result)

    print("Finished processing individual tables")
    df_list = [df for df in results if df is not None]

    if df_list:
        print('Concatanating dataframes' )
        all_df = pd.concat(df_list, axis=0, ignore_index=True)
        print("Finished processing LAST tables")
        return all_df
    else:
        print("No valid data found in LAST files.")
        return pd.DataFrame()

def lenset(list_):
    return len(set(list_))

def replace(data, level, name, new_name):
    group = set(data[data[level] == name]['Family'])
    return {x: new_name for x in group}

def replacement_process(df_top, replacement_df):
    print("Start replacing process...")
    list_replacement = [
        replace(df_top, row['level'], row['old_name'], row['new_name'])
        for _, row in replacement_df.iterrows()
    ]

    super_dict = {key: val for d in list_replacement for key, val in d.items()}
    large_list = set(super_dict.keys())
    rest_species = {x: "Other eukaryotes" for x in df_top['Family'].unique() if x not in large_list}

    super_dict.update(rest_species)
    return super_dict

def dedup_df(df):
    df[["subject_id", "rrna_type", "rrna_loc", "clone", "Domain","Clade", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"]] = df['subject id'].str.split("|", expand=True)
    df = df[df['rrna_type'] == "18S_rRNA"].copy()
    
    
    # Get only the best hit for each read
    df['rank'] = df.groupby(['query id'])['bit score'].rank("min", ascending=False)
    df = df[df['rank'] == 1]
    return df

def filter_df(df):
    global identity, length
    df = df[df['% identity'] >= identity]
    df = df[df['alignment length'] >= length]
    print(lenset(df['cell']), "cells passed identity and alignment length threshold")
    return df

def edit_df(filtered_df):
    global sum_df
    sum_df = filtered_df.groupby(['cell', level])['query id'].agg('count').reset_index()
    sum_df = sum_df[sum_df['query id'] >= count]
    sum_df['cell'] = sum_df['cell'].apply(lambda x: '.'.join(x.split('.', 2)[:2]))
    sum_df['percent'] = sum_df['query id'] / sum_df.groupby('cell')['query id'].transform('sum')
    return sum_df

def load_last_tables(last_path):
    last_df = load_last(last_path)
   # last_df.to_csv('18s_last_df.csv', index=False) # We don't need it since it's faster
    return last_df

def process_file(file, last_path):
    """Process a single LAST file and return the deduplicated dataframe."""
    parts = file.split('.', 2)
    cell = '.'.join(parts[:2]).replace(last_path, "")
    
    columns = ["query id", "subject id", "% identity", "alignment length", 
               "mismatches", "gap opens", "q. start", "q. end", 
               "s. start", "s. end", "evalue", "bit score"]

    last_df = pd.read_csv(file, skip_blank_lines=False, sep='\t', header=None, names=columns)
    
    last_df = last_df[~last_df['query id'].str.contains(r'#', na=True)]
    last_df['cell'] = cell

    if lenset(last_df['query id']) > 0:
        return dedup_df(last_df)
    return None  # Skip files with no valid data

def get_args():
    parser = ArgumentParser(description="Choose a subset of cells that pass a criteria")
    parser.add_argument("--last_path", required=True, help="Path for all the LAST files")
    parser.add_argument("--last_summary", default='False', help="Path for the summarized LAST table")
    parser.add_argument("--replacement_file", required=True, help="A tab-delimited file for replacing taxonomic names")
    parser.add_argument("--identity", default=99, type=float, help="Threshold for read percent homology")
    parser.add_argument("--length", default=85, type=float, help="Threshold for read alignment length")
    parser.add_argument("--count", default=10, type=float, help="Threshold for number of reads to consider")
    parser.add_argument("--level", default='Annotation', help="Taxonomic level for data subset")
    parser.add_argument("--workers", default=4, type=int, help="Max number of parallel processes")
    return parser.parse_args()

def main():
    args = get_args()
    
    global identity, length, count, level, max_workers
    identity, length, count, level, max_workers = args.identity, args.length, args.count, args.level, args.workers
    
    last_path = args.last_path
    last_summary = args.last_summary
    replacement_file = args.replacement_file

    if last_summary.lower() != 'false':
        last_df = pd.read_csv(last_summary)
    else:
        last_df = load_last_tables(last_path)

    replacement_df = pd.read_table(replacement_file, sep="\t")
    super_dict = replacement_process(last_df, replacement_df)

    last_df['Annotation'] = last_df['Family'].map(super_dict)
    df_dedup = last_df.drop_duplicates(subset=["query id", level]).drop_duplicates('query id', keep=False)
    edited_df = filter_df(df_dedup)

    final_df = edit_df(edited_df)
    final_df.to_csv(f'last_summary.{level}.csv', index=False)

if __name__ == "__main__":
    main()
