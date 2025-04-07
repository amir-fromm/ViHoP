#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
from tqdm import tqdm
from argparse import ArgumentParser


# In[8]:


def edit_blast_table(tmp_df, dict_annotations):
    tmp_df['cell'] = tmp_df['qseqid'].str.split(r'i\d_',expand = True)[1]
    tmp_df['Annotation'] = tmp_df['cell'].map(dict_annotations)
    
    tmp_df['Annotation'] = tmp_df['Annotation'].fillna('Unknown')
    return tmp_df

def subset_blast_table(tmp_df, identity, length):
    clean_df = tmp_df.copy()

    print("Filtering by identity and alignment length...")
    clean_df = clean_df[(clean_df['pident'] >= identity) & (clean_df['length'] >= length)]
    print("Choosing only best hit...")
    clean_df['rank'] = clean_df.groupby(['qseqid'])['bitscore'].rank("min",ascending = False)
    clean_df = clean_df[clean_df['rank'] == 1]
   # clean_df = clean_df[clean_df['Annotation'] != 'Metazoa']
    return clean_df


# In[9]:


def load_data(annotation_file,blast_file,taxa_file):
    annotation = pd.read_csv(annotation_file)
    
    list_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(" ")
    blast = pd.read_table(blast_file, sep = '\t', names = list_columns)
    
    taxa = pd.read_csv(taxa_file, header = None, names = ['Virus', 'Family'])
    return annotation, blast, taxa


# In[10]:


def virus_taxonomy(new_df, dict_taxonomy):
    with tqdm(total=new_df.shape[0]) as pbar:
        for index, row in new_df.iterrows():
            subject = row['sseqid']
          #  new_subject = subject.replace("-", "_")
            matched = False  # Flag to check if any key matched
            for key in dict_taxonomy.keys():
                if key in subject or subject in key:
                    new_df.loc[index, 'Virus family'] = dict_taxonomy[key]
                    matched = True  # Set the flag when a match is found
                    break  # Optional: break if you only need one match, remove if multiple matches are handled differently

            if not matched:  # If no key matched after checking all keys
                new_df.loc[index, 'Virus family'] = 'incertae_sedis'
            pbar.update(1)
    return new_df


# In[11]:


# Saving files
def saving_files(output_path,output_name,tmp_df):
    with open (output_path + "/virus_contigs.txt", "w") as f:
        for seq in tmp_df['qseqid'].unique():
            f.write(f"{seq}\n")
            
    tmp_df.to_csv(output_path + "/" + output_name, sep = '\t', index = False)


# In[ ]:


def get_args():
    parser = ArgumentParser(description="Choose a subset of cells that pass a criteria")
    parser.add_argument("--output_path", help = "Path for the output files", required=True)
    parser.add_argument("--taxa_file", help = "Path for taxonomy file", required=True)
    parser.add_argument("--blast_file", help = "Path for all the virus Blast files", required=True)
    parser.add_argument("--annotation_file", help = "Path for all the last files", required=True)
    parser.add_argument("--identity", help = "Threshold for read percent homology", default=90, type = float)
    parser.add_argument("--length", help = "Threshold for alignment length", default=100, type = float)
    parser.add_argument("--output_name", help="Name of output file", default = 'summary_blast_viruses.tsv')

    return parser.parse_args() 


# In[18]:


def main():
    args = get_args()
    print('Summarizing virus blast results...')
    annotation_file = args.annotation_file
    blast_file = args.blast_file
    taxa_file = args.taxa_file
    output_path = args.output_path
    output_name = args.output_name
    identity = args.identity
    length = args.length
    
    annotation,blast,taxa = load_data(annotation_file,blast_file,taxa_file)

    dict_taxonomy = dict(zip(taxa['Virus'],taxa['Family'])) 
    dict_annotations = dict(zip(annotation['cell'],annotation['Blast Annotation']))


    blast = edit_blast_table(blast,dict_annotations)
    
    blast_df = subset_blast_table(blast, identity,length)
    blast_df['Virus family'] = ''
    new_df = blast_df.copy(deep = True)

    summary_df = virus_taxonomy(new_df, dict_taxonomy)

    saving_files(output_path,output_name,summary_df)

if __name__ == "__main__":
    main()