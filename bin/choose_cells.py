#!/usr/bin/env python
# coding: utf-8

# In[46]:


from argparse import ArgumentParser



# In[47]:


import scanpy as sc
import os
import pandas as pd


# In[48]:


def load_cellranger_matrix(path, sample_name=None):
    """
    Load a 10x Genomics Cell Ranger matrix into an AnnData object.

    Parameters:
    - path (str): Path to directory containing 'matrix.mtx.gz', 'features.tsv.gz', and 'barcodes.tsv.gz'
                  Usually this is the 'filtered_feature_bc_matrix' or 'raw_feature_bc_matrix' directory.
    - sample_name (str): Optional. If provided, assigns sample name to `.obs['sample']`.

    Returns:
    - AnnData object with raw counts in .X
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Path does not exist: {path}")
    print(path)
    # Load using Scanpy’s built-in reader
    adata = sc.read_10x_mtx(path, var_names='gene_symbols')

    # Assign sample metadata if provided
    if sample_name:
        adata.obs["sample"] = sample_name

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    return adata


# In[49]:


def concatenate(sample_adata_list,sample_names):
    adata = sample_adata_list[0].concatenate(
        sample_adata_list[1:], 
        batch_key="sample",
            join="outer",
        batch_categories=sample_names,
        index_unique=None  # ← keep original obs_names from each sample
    )
    return adata


# In[50]:


def get_all_cells_from_anndata(adata, formatted = False):
    # Get all cells from an unforamtted anndata (without formatted cells names) or from a formatted anndata
    if formatted == True:
        list_cells = list(adata.obs.index)
    else:
        list_cells = list(adata.obs.reset_index()['index'].str.split("-1", expand = True)[0] + "." + adata.obs.reset_index()['sample'])
    return list_cells


# In[51]:


def load_data(paths):
  

    adatas = [load_cellranger_matrix(p, sample_name=s) for s, p in paths.items()]

    sample_names = paths.keys()

    adata = concatenate(adatas,sample_names)
    adata.obs_names = [
        f"{bc.replace('-1','')}.{samp}" for bc, samp in zip(adata.obs_names, adata.obs["sample"])
    ]
    return adata


# In[52]:


def get_args():
    parser = ArgumentParser(description="Parallel BLAST summarizer for large files")
    parser.add_argument("--sample_list", required=True, help="List of samples")
    parser.add_argument("--path", required=True, help="Path to sample and output")
    parser.add_argument("--sample_table", required=True, help="A table of samples and corresponding Fastq files")
    parser.add_argument("--min_gene", required=False, default = 1, help="Minimum cells expressing a gene (Sum of cells expressing genes). default=1")
    parser.add_argument("--min_umi", required=False, default = 1, help="Minimum expression per cell (sum of UMIs). default=1")
    parser.add_argument("--output", required=False, default = 'list_cells', help="Output file name. default=list_cells.tsv")

    return parser.parse_args()


# In[53]:


def main():
    args = get_args()
    sample_list = args.sample_list.split(",")
    path = args.path
    paths = {sample:os.path.join(path,sample,'outs/filtered_feature_bc_matrix') for sample in sample_list}
    adata = load_data(paths)
    min_gene = args.min_gene
    min_umi = args.min_umi
    print("Minimum cells expressing gene:",min_gene,"\nMinimum UMI per cell:", min_umi)
    sc.pp.filter_genes(adata, min_counts=int(args.min_gene))
    sc.pp.filter_cells(adata, min_counts=int(args.min_umi))
    print(adata)
    total_cells = get_all_cells_from_anndata(adata,formatted = True )
    table_samples = pd.read_table(args.sample_table,sep = "\t", header = None,index_col = 0)
    table_samples = table_samples.squeeze()
    dict_samples = table_samples.to_dict()
    df = pd.DataFrame([v.split('.') for v in total_cells])
    df.columns = ["cell", "sample"]
    df['fastq'] = df['sample'].map(dict_samples)
    output_file = args.output
    df[['cell','fastq']].to_csv(os.path.join(path, output_file), sep = "\t",index = False,header=None)


# In[54]:


if __name__ == "__main__":
    main()

