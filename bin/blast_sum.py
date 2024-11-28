#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(description="Choose a subset of cells that pass a criteria")
    parser.add_argument("--path", help = "Path for all the files", required=True)

    parser.add_argument("--metapr2_file",help = "metapr2_file = The output files for metapr2 blast file", required = True)
    parser.add_argument("--pr2_file",help = "pr2_file =  The output files for pr2 blast file" , required=True)
    parser.add_argument("--replacement_file",help = "A tab delimited file for replacing taxonomic names and assigning colors for each taxonomic group", required=True)
  
    return parser.parse_args()


# Dedup functions remove conflicting contigs (contigs with the best match for multiple annotations)
# and drop ambiguous cells (cells with contigs matching 18S rRNA of multiple organisms)
def dedup_data(path,file,data_type,replacement_df,cells):
    data = load_data(path,file,data_type)
    global super_dict
    super_dict = replacement_process(data,replacement_df)
    data['Annotation'] = data['Family'].map(super_dict)
    data['barcode'] = data['cell'].str.split(".", expand = True)[0]
    if cells != -1:
        data = data[data['cell'].isin(cells)]
    data = subset_results(data)
    data = data.drop_duplicates(subset = ["qseqid",'Annotation']).drop_duplicates('qseqid',keep = False)
    data = get_chimeric_cells(data)

    return data

def load_data(path,file,db):
    list_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(" ")

    print("reading file " + file + "...")
    data = pd.read_table(path+file, sep = "\t", names = list_columns)
    data['cell'] =  data['qseqid'].str.split("i\d+_", expand = True)[1]
    # in case of pr2 open:
    if db == "PR2":
        data[["sseqid","rrna_type","rrna_loc","clone","Domain","Clade","Supergroup","Division","Class","Order","Family","Genus","Species"]] = data['sseqid'].str.split("|",expand = True)
        data = data[(data['rrna_type'] == "18S_rRNA")]
  
    # in case of metapr2 open:
    elif db == "metaPR2":
        data[["sseqid","Clade","Supergroup","Division","Class","Order","Family","Genus","Species"]] = data['sseqid'].str.split("|",expand = True)
    else:
        print("Unrecognized data type! Exiting...")
        exit()
    print("Finished processing file")
    return data


# In[25]:


def subset_results(sample_data):
    
    print("Filtering queries by by length...")
    sample_data = sample_data[(sample_data['length'] >= 100)]
    
    print("sorting data...")
    sample_data = sample_data.sort_values(by = ['qseqid','bitscore'],ascending = False)
    
    print("Choosing only best hit...")
    sample_data['rank'] = sample_data.groupby(['qseqid'])['bitscore'].rank("min",ascending = False)
    sample_data = sample_data[sample_data['rank'] == 1]
    
    print("Filtering queries by identity...")
    sample_data = sample_data[(sample_data['pident'] >= 99)]
    return sample_data


def get_chimeric_cells(df):
    print("Find chimeric cells...")
    chimeric_contigs = find_chimeras(df)
    chimeric_cells = set(df[df['qseqid'].isin(chimeric_contigs)]['cell'])
    print('Found', len(chimeric_cells), 'chimeric cells')
    return df[~df['cell'].isin(chimeric_cells)]

def lenset(x):
    return len(set(x))




#create annotation data
def annotation_data(pr2_data,metapr2_data):
    print("Annotating cells")
    global data_annotations
    # create the data annotations dataset
    sub_pr2 = pr2_data[['qseqid','cell','Annotation','length']]
    sub_pr2['cell'] = sub_pr2['cell'].str.strip()
    
    sub_metapr2 = metapr2_data[['qseqid','cell','Annotation','length']]
    sub_metapr2['cell'] = sub_metapr2['cell'].str.strip()
    
    data_annotations = pd.concat([sub_pr2,sub_metapr2])
    
    
    # remove ambigous contigs between the databases
    data_annotations = data_annotations.sort_values(by = 'length', ascending = False)
        
    data_annotations = data_annotations.drop_duplicates(subset = ["qseqid",'Annotation'])
    

    data_annotations['percent'] = data_annotations['length'] / data_annotations.groupby('cell')['length'].transform('sum')
    data_annotations = data_annotations.groupby(['cell','Annotation',]).agg({'percent':'sum','length':'max'}).reset_index()
    return data_annotations


# Create a dictionary for each of the taxonomic groups listed in the input file
def replace(data,level, name, new_name):
    group = set(data[data[level] == name]['Family'])
    value_to_replace = {x : new_name for x in group}
    return value_to_replace

# Iterate through all the main taxonomic groups listed in the input file and unify each of them under the name listed as replacement
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
    print ("Finished replacing process")
    return super_dict


def find_chimeras(df):
    df = df.drop_duplicates(subset = ["qseqid",'Annotation'])
    # Sort the dataframe by 'Gene ID' and 'Col1' for easier processing
    df = df.sort_values(by=['qseqid', 'qstart'])

    # Initialize lists to store the results
    gene_ids = []
    col1_values = []
    col2_values = []
    collapsed_ranges = []
    individual_range_counts = []

    current_gene_id = None
    current_col1 = None
    current_col2 = None
    current_range = None
    current_individual_ranges = set()

    for index, row in df.iterrows():
        gene_id = row['qseqid']
        col1 = row['qstart']
        col2 = row['qend']

        if current_gene_id is None:
            current_gene_id = gene_id
            current_col1 = col1
            current_col2 = col2
            current_range = f"{col1}-{col2}"
            current_individual_ranges.add(current_range)
        elif gene_id == current_gene_id and col1 <= current_col2:
            current_col2 = max(current_col2, col2)
            current_range = f"{current_col1}-{current_col2}"
            current_individual_ranges.add(f"{col1}-{col2}")
        else:
            gene_ids.append(current_gene_id)
            col1_values.append(current_col1)
            col2_values.append(current_col2)
            collapsed_ranges.append(current_range)
            individual_range_counts.append(len(current_individual_ranges))

            current_gene_id = gene_id
            current_col1 = col1
            current_col2 = col2
            current_range = f"{col1}-{col2}"
            current_individual_ranges = {current_range}

    # Append the last group
    if current_gene_id is not None:
        gene_ids.append(current_gene_id)
        col1_values.append(current_col1)
        col2_values.append(current_col2)
        collapsed_ranges.append(current_range)
        individual_range_counts.append(len(current_individual_ranges))

    # Create a new dataframe with the collapsed ranges and individual range counts
    collapsed_df = pd.DataFrame({'qseqid': gene_ids, 'qstart': col1_values, 'qend': col2_values, 'Collapsed Ranges': collapsed_ranges, 'Individual Range Counts': individual_range_counts})

    
    collapsed_df_new = collapsed_df[collapsed_df['Individual Range Counts'] >= 1]
    if len(collapsed_df_new) > 0:
        dict_ranges_count = dict(collapsed_df_new.groupby('qseqid')['Collapsed Ranges'].apply(list).str.len())

        collapsed_df_new['Ranges count'] = collapsed_df_new['qseqid'].map(dict_ranges_count)
        chimeric_contig =  set(collapsed_df_new[(collapsed_df_new['Ranges count'] >= 2)]['qseqid'])    
    else:
        chimeric_contig = {}
    return chimeric_contig


def main(path,
         metapr2_file,
         pr2_file,
        replacement_file):
 
    cells = -1
    names_df = pd.read_table(path+ "/" + replacement_file, sep = "\t")
    if 'metapr2_data' not in globals() or 'pr2_data' not in globals():
        global metapr2_data, pr2_data
        metapr2_data = dedup_data(path = path,
            file = metapr2_file , data_type = "metaPR2",replacement_df = names_df, cells = cells)
        pr2_data = dedup_data(path = path,
            file = pr2_file , data_type = "PR2",replacement_df = names_df,cells = cells)

    data_annotations = annotation_data(pr2_data,metapr2_data)

    data_annotations.to_csv(path + "/annotation_permissive.csv")
    data_annotations.drop_duplicates(['cell','Annotation']).groupby('cell')['Annotation'].apply(list).to_csv(path + "/annotation_permissive_groups.csv")



# In[32]:


'''
Input:
# path = Path for all the files
# metapr2_file = The summary output files for metapr2 blast file (output of the blast summary script)
# pr2_file =  The summary output files for pr2 blast file (output of the blast summary script)
# replacement_file = A tab delimited file for replacing taxonomic names and assigning colors for each taxonomic group.
Must contain most of the main groups present in the analysis (an inspection of the pr2 and metapr2 results is required)
Contains the following columns:
    old_name: The taxonomic name to replace (For example: Bacillariophyta)
    level: The taxonomic level of the old name (For example: Class)
    new_name: The new name to use in the plot (For example: Diatoms)
    color: The color assigned to the taxonomic group in the Sankey plot (For example: #8C613C)
# output_file = The name of the output png image (For exmaple: sankey.png)

'''
if __name__ == "__main__":
   args = get_args()
   main(args.path,
        args.metapr2_file,
        args.pr2_file,
        args.replacement_file)




