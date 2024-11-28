#!/usr/bin/env python
# coding: utf-8

# In[14]:


import glob
import pandas as pd
from tqdm import tqdm
from argparse import ArgumentParser



# Load all last output tables, combine them and get only best hit for each read
def load_last(last_path):
    print("Processing last tables")
    if not 'all_df' in globals():
        global df_list
        df_list = []
        no_list = []
        flist = glob.glob(last_path + "*.last")
        with tqdm(total=len(flist)) as pbar:
            for file in flist:
                parts = file.split('.', 2)
                # Join the first two parts with a dot
                cell = '.'.join(parts[:2])
                cell = cell.replace(last_path,"")
                columns = "query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score".split(", ")
                last_df = pd.read_csv(file, skip_blank_lines = False, sep = '\t',header=None, names = columns)

                last_df = last_df[~last_df['query id'].str.contains(r'#', na=True)]

                last_df['cell'] = cell

                if lenset(last_df['query id']) > 0:
                    short_df = dedup_df(last_df)


                    df_list.append(short_df)
                else:
                    no_list.append(last_df)
                pbar.update(1)
        all_df = pd.concat(df_list, axis=0, ignore_index=True)
    
    
    print("Finished processing last tables")
    return all_df



def lenset(list_):
    return len(set(list_))


def replace(data,level, name, new_name):
    group = set(data[data[level] == name]['Family'])
    value_to_replace = {x : new_name for x in group}
    return value_to_replace

def replacement_process(df_top,replacement_df):
    print ("Start replacing process...")
    list_replacement = []
    for index, row in replacement_df.iterrows():
        list_replacement.append(replace(df_top,row['level'], row['old_name'],row['new_name']))
      
    super_dict = {key:val for d in list_replacement for key,val in d.items()}
    large_list = set(super_dict.keys())
 
    rest_species = {x : "Other eukaryotes" for x in [x for x in list(df_top.Family) if x not in large_list]}
    
    dicts = [super_dict,
    rest_species]
    super_dict = {key:val for d in dicts for key,val in d.items()}
    return super_dict


def dedup_df(df):
   # print(lenset(df['cell']), "cells is total")
    
    filtered_df_copy = df.copy()
    # Split the 'subject id' column and assign the results to new columns
    filtered_df_copy[["subject_id", "rrna_type", "rrna_loc", "clone", "Domain","Clade", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"]] = filtered_df_copy['subject id'].str.split("|", expand=True)
    filtered_df_copy = filtered_df_copy[(filtered_df_copy['rrna_type'] == "18S_rRNA")]
    
    # Get only the best hit for each read
    filtered_df_copy['rank'] = filtered_df_copy.groupby(['query id'])['bit score'].rank("min",ascending = False)
    filtered_df_copy = filtered_df_copy[filtered_df_copy['rank'] == 1]
    
    
 
    return filtered_df_copy




def filter_df(df, identity, length):
# Remove reads that are mapped to different classes
    
    
    df = df[df['% identity'] >= identity]
    df = df[df['alignment length'] >= length]
    
    print(lenset(df['cell']), "cells passed identity and alignment length threshold")
    return df


# In[21]:


def edit_df(filtered_df):
    global sum_df
    sum_df = filtered_df.groupby(['cell','Annotation'])['query id'].agg('count').reset_index()
    

    sum_df['cell'] = sum_df['cell'].apply(lambda x: '.'.join(x.split('.', 2)[:2]))

    sum_df['percent'] = sum_df['query id'] / sum_df.groupby('cell')['query id'].transform('sum')


    
    return sum_df


# In[22]:


def load_last_tables(last_path, output_path):
    
    last_df = load_last(last_path)
    last_df.to_csv(output_path + '/18s_last_df.csv')
    return last_df



def get_args():
    parser = ArgumentParser(description="Choose a subset of cells that pass a criteria")
    parser.add_argument("--output_path", help = "Path for the output files", required=True)
    parser.add_argument("--last_path", help = "Path for all the last files", required=True)
    parser.add_argument("--last_summary", help = "Path for the summarized last table", default='False')
    parser.add_argument("--replacement_file",help = "A tab delimited file for replacing taxonomic names and assigning colors for each taxonomic group", required=True)
    parser.add_argument("--identity", help = "Threshold for read percent homology", default=99, type = float)
    parser.add_argument("--length", help = "Threshold for read alignment length", default=85, type = float)

    return parser.parse_args()

def main(last_path,
last_summary,
replacement_file,
output_path):
    

    if last_summary != 'false':

        last_df = pd.read_csv(last_summary)

    else:
        last_df = load_last_tables(last_path, output_path)


    replacement_df = pd.read_table(replacement_file, sep = "\t")

    super_dict = replacement_process(last_df,replacement_df)

    last_df['Annotation'] = last_df['Family'].map(super_dict)

    df_dedup = last_df.drop_duplicates(subset = ["query id",'Annotation']).drop_duplicates('query id',keep = False)
    edited_df = filter_df(df_dedup, identity, length)

    final_df = edit_df(edited_df)

    final_df.to_csv(output_path + '/18s_mapping_sum.csv')


if __name__ == "__main__":
    args = get_args()
    last_path = args.last_path
    last_summary = args.last_summary
    replacement_file = args.replacement_file
    output_path = args.output_path
    identity = args.identity
    length = args.length
    main(last_path,
         last_summary,
         replacement_file,
         output_path)
