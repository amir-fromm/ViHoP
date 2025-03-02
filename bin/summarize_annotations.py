#!/usr/bin/env python
# coding: utf-8

# Load the list of infected cells, the result of the blast search and the result of the 18s last search.
# Returns a combined dataframe of the results, and a summarized dataframe of the cells the pass all thresholds and not chimeric

import pandas as pd
from argparse import ArgumentParser

def load_data(path_cells, path_last,path_blast):
    list_cells = pd.read_table(path_cells, names = ['cell', 'path'])
    last_df = pd.read_csv(path_last, index_col=0)
    blast_df = pd.read_csv(path_blast, index_col=0)
    return list_cells,last_df,blast_df

# Create the metadata columns for the list of infected cells
def edit_list_cell(list_cells):
    list_cells['Sample'] = list_cells['path'].str.split('/').str[-1]
    list_cells['Barcode'] = list_cells['cell']
    list_cells['cell'] = list_cells['Barcode'] +"."+ list_cells['Sample'] 
    return list_cells

# Edit the dataframes so they will have similar shared columns and different names for columns that shared the same name
def edit_dfs(list_cells,last_df,blast_df):
    
    df_cells =  edit_list_cell(list_cells)
    blast_df[['Barcode','Sample']] = blast_df['cell'].str.split('.',expand = True)
    last_df[['Barcode','Sample']] = last_df['cell'].str.split('.',expand = True)
    last_df = last_df.rename(columns = {'Annotation':'last Annotation','query id':'last # mapped reads','percent':'last percent mapped reads'})
    blast_df = blast_df.rename(columns = {'Annotation':'Blast Annotation','length':'Blast max alignment length','percent':'Blast percent alignment' })
  
    return list_cells,last_df,blast_df

# merges the dataframes.
# Filter results by creteria: at least 10 mapped cells (last), at least 90% of reads match one taxa (last), at least 90% of alignment length match one taxa (blast)
# checks if last annotation equals blast annotation
def summarize_df(df_cells,last_df,blast_df,count,percent):
    merge_on = ['cell','Sample','Barcode']
    full_db = pd.merge(pd.merge(blast_df, last_df, on = merge_on,how='outer'),df_cells[['cell','Sample','Barcode']] ,on = merge_on,how='outer')
    # Check if all pass threshold
    sum_df = full_db[(full_db['last # mapped reads'] >= count )  & (full_db['last percent mapped reads'] >= percent) & (full_db['Blast percent alignment'] >= percent)]
    sum_df = sum_df[sum_df['Blast Annotation'] == sum_df['last Annotation']]
    return full_db,sum_df

def write_files(full_db, sum_df, output_path):
    full_db.to_csv(output_path + '/full_annotation.csv')
    sum_df.to_csv(output_path + '/summarized_annotation.csv')

def get_args():
    parser = ArgumentParser(description="Choose a subset of cells that pass a criteria")
    parser.add_argument("--output_path", help = "Path for the output files", required=True)
    parser.add_argument("--last_path", help = "Path for all the last files", required=True)
    parser.add_argument("--blast_path", help = "Path for all the Blast files", required=True)
    parser.add_argument("--cell_path", help = "Path for all the last files", required=True)
    parser.add_argument("--count_reads", help = "Threshold for minimum number of mapped reads to have annotation", required=True) 
    parser.add_argument("--percent", help = "Threshold for percent mapped reads/alignment length to determine if chimeric or not", required=True) 
    return parser.parse_args() 

def main():
    args = get_args()
    print('Summarizing annotation results...')
    list_cells,last_df,blast_df = load_data(args.cell_path,
                                            args.last_path,
                                            args.blast_path)
    print('Editing annotation tables')
    list_cells,last_df,blast_df = edit_dfs(list_cells,last_df,blast_df)
    print('Summerizing annotation tables')
    count=int(args.count_reads)
    percent=int(args.percent)
    percent=percent/100
    full_db, sum_df = summarize_df(list_cells,last_df,blast_df,count,percent)
    print('Writing results:', sum_df.shape[0], 'cells were annotated and are not chimeric')
    write_files(full_db, sum_df, args.output_path)
    print('Done')
if __name__ == "__main__":
     main()