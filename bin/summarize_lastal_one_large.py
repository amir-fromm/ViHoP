#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from argparse import ArgumentParser

import pandas as pd
import subprocess, shlex
import subprocess
from tqdm import tqdm



# In[ ]:


# get the nodes data from the refseq db
# make a dictionary from the taxonomy ID and the division


def load_merged():
    path = '/home/labs/vardi/Collaboration/databases/refseq/'
    merged_file = 'merged.edited.dmp'
    merged_df = pd.read_csv(path + merged_file,sep = '\t',header = None,names = ['old','new'],index_col = False)
    dict_merged = dict(zip(merged_df['old'].astype(str),merged_df['new'].astype(str)))
    return dict_merged
# In[ ]:


# Load all last output tables from one large table in chunks, get only best hit for each read and add to a list
def load_last(last_path, chunksize, total_lines):
    print("Processing large last table")
    dtype_mapping = {
    "query id": "str",
    "subject id": "str",
    "% identity": "float64",
    "alignment length": "int64",
    "bit score": "float64",  # ✅ Explicitly set float64 to avoid dtype issues
    }
    global df_list
    df_list = []
    columns = "query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score".split(", ")
    
        
    for chunk in tqdm(pd.read_csv(last_path, skip_blank_lines = False, comment="#",  sep = '\t',header=None, names = columns, chunksize=chunksize), 
                  total=total_lines // chunksize + 1, desc="Processing Chunks"):
        
        chunk["rank"] = chunk.groupby("query id")["bit score"].rank("min", ascending=False)
        
        #idx = chunk.groupby("query id")["bit score"].idxmax()
      #  best_hits_chunk = chunk.loc[idx]
        
        df_list.append(chunk[chunk["rank"] == 1])
    last_df = pd.concat(df_list)
    print("Ranking again to to avoid conflicting chunks")
    last_df["rank"] = last_df.groupby("query id")["bit score"].rank("min", ascending=False)
    last_df = last_df[last_df["rank"] == 1]
  #  last_df = all_df.sort_values('bit score', ascending = False).drop_duplicates('query id', keep = False)
    print("Finished processing last table")
    return last_df


# Define column names
columns = ["query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens",
           "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"]

chunksize = 1_000_000  # Adjust based on RAM
df_list = []


# write the hit files as a temporary file
def write_hits(last_df):
    hits = set(last_df['subject id'])
    output = '/home/labs/vardi/amirz/tmp/tmp_hits.txt'
    with open(output,"w") as f:
        for hit in hits:
            f.write(hit + '\n')
    return output


# In[ ]:


# Using grep to get the taxonomy ID from all the hits from the temporary file. 
# Turn the list into a dictionary of protein and taxonomy ID
def hit_table(output):
    catalog_file = '/home/labs/vardi/Collaboration/databases/refseq/RefSeq-release223.taxid.catalog'

    cmd = "grep -f /home/labs/vardi/amirz/tmp/tmp_hits.txt " + catalog_file
    cmd2 = shlex.split(cmd)
    proc = subprocess.Popen(cmd2, stderr = open("error3.txt", "w"), stdout=subprocess.PIPE)
    
    output = proc.stdout.read().decode('utf-8')
    hits_list = output.split("\n")
   # print(hits_list)
    table = [x.split(" ") for x in hits_list if len(x) > 1]
    tmp_df = pd.DataFrame(table)
    dict_proteins = dict(zip(tmp_df[1],tmp_df[0]))

    return dict_proteins
#end = time.time()



# In[ ]:


# Assign domain (Eukryoters/Bacteria/Viruses/Phages/Others) to TaxID through the divisions
def assign_division(dict_proteins,dict_nodes,dict_domains,dict_merged):
    print("Assigning division to proteins")
    protein_to_division = {}
    for protein, taxid in dict_proteins.items():
        try:
            protein_to_division[protein] = dict_domains[dict_nodes[taxid]]
        except:
            if taxid in dict_merged.keys():
                new_taxid = dict_merged[taxid]
                print ("Merged taxid", taxid, ":", new_taxid)
                if new_taxid in dict_nodes.keys():
                    protein_to_division[protein] = dict_domains[dict_nodes[new_taxid]]
                else:
                    print(new_taxid, "taxid for protein",protein, "not in nodes dictionary")
            else:
                print ("Unknown taxid", taxid, "for protein", protein)
                protein_to_division[protein] = "Unknown taxid"
    return protein_to_division


# In[ ]:


# Edit the output dataframe by assigning taxonomy to the viruses using Aylward DB
def edit_last_df(last_df, protein_to_division):
    print("Editing the output dataframe")
    last_df['Domain'] = last_df['subject id'].map(protein_to_division)
    last_df['Domain'] = last_df['Domain'].fillna('Viruses')
    

    # Drop conflicting domains
    final_df = last_df.drop_duplicates(subset = ["query id",'Domain']).drop_duplicates('query id',keep = False)
    
    return final_df

# In[ ]:


def load_nodes():
    path = '/home/labs/vardi/Collaboration/databases/refseq/'
    nodes_file = 'nodes.short.dmp'
    global nodes_df
    nodes_df = pd.read_csv(path + nodes_file,sep = ' ',header = None,names = ['taxid','Division'],index_col = False)
    
    dict_nodes = dict(zip(nodes_df['taxid'].astype(str),nodes_df['Division']))
    return dict_nodes

def save_files(output_path,new_df):
    print("Saving files")
    new_df.to_csv(output_path , sep = '\t', index=False)

def main():
    args = get_args()
    last_path = args.last_path
    chunks = args.chunks
    size_table = args.size_table
    output_path = args.output_path    
    chunks = args.chunks
    size_table = args.size_table

    if 'last_df' not in globals():
        global last_df

        last_df = load_last(last_path, 
                            chunks, 
                            size_table)
    
        
    output = write_hits(last_df)
    dict_proteins = hit_table(output)
    #list_taxids = set(tmp_df[0])
    dict_nodes = load_nodes()
    dict_merged = load_merged()
   ### dict_virus_class = virus_tax_data()
    dict_domains = {0: 'Bacteria',
     1: 'Eukaryotes',
     2: 'Eukaryotes',
     3: 'Phages',
     4: 'Eukaryotes',
     5: 'Eukaryotes',
     6: 'Eukaryotes',
     7: 'Others',
     8: 'Others',
     9: 'Viruses',
     10: 'Eukaryotes',
     11: 'Others'}
    
    protein_to_division = assign_division(dict_proteins,dict_nodes,dict_domains,dict_merged)
    global new_df
    new_df = edit_last_df(last_df, protein_to_division)

    save_files(output_path,new_df)
    print("Done!")


def get_args():
    parser = ArgumentParser(description="Choose a subset of cells that pass a criteria")

    
    parser.add_argument("--last_path", help = "Path for the output lastal file", required=True)
    parser.add_argument("--chunks", help = "Size of chunks to take from table. default = 10000", default=10000, type = float)
    parser.add_argument("--size_table", help = "The total number of rows in the last table", required=True, type = float)
    parser.add_argument("--output_path", help="Absolute path of llastalst results", default = 'summary_last.tsv')

    return parser.parse_args() 

if __name__ == '__main__':
    main()