  
# ViHoP: Virus-Host Pairing
## Pairing viruses to their hosts from environmental single-cell 10x data
Code repository and instructions for finding host-virus pairs from single-cell 10X sequencing data

amir.fromm@weizmann.ac.il  
## What the script does
This pipeline takes single-cell RNA-sequenceing data 10X sequecing and using a viral and host references, produce a matrix of host cells and the viruses that infects them.
This repository is built to work on an LSF platform since many of the script use parallization.

![alt text](https://github.com/amir-fromm/ViHoP/blob/main/sc_virus_host_detection_pipeline.png?raw=true)

### 1. Cloning this repository

```
git clone https://git
``` 
Then go into the directory
```
cd host-virus pairing
```
Create environment with the requirement file
```
conda create --name host-virus-pairing --file requirements.txt
conda activate host-virus-pairing
```

### 2. Other dependencies  

Besides the packages from the requirements file, the pipeline also uses modules that were installed in our lsf platform and are called module load command:
```
CD-HIT/4.6.6
CellRanger/8.0.0
TrimGalore/0.6.5
SPAdes/3.15
BBMap
BLAST/2.11.0
seqtk

```
The speicfic module names can be changed using the config file. If you don't have these tools installed on your lsf platform, install them on your conda environment.

### 3. Making host reference databases

You need to to download the latest protist host 18s rRNA databases PR2 (version 5) and metaPR2 (version 2) (https://pr2-database.org/, https://app.metapr2.org/metapr2/).
You need to download these databases and create last and blast db from them:
```
makeblastdb -in db_path/PR2.fasta -blastdb_version 5 -title "PR2" -dbtype nucl -out PR2
makeblastdb -in db_path/metaPR2.fasta -blastdb_version 5 -title "metaPR2" -dbtype nucl -out metaPR2
lastdb -c -P 4 PR2 db_path/PR2.fasta

```
### 4. Virus taxonomy

You need to to have a virus database of your choice (you can use Aylward's giant virus database: https://faylward.github.io/GVDB/), and a virus taxonomy table should be available. 
It is a comma-delimited file where the first column is the contig name in the virus reference (or part of it) and the second column is the virus taxonomy.
Afterv you download the virus database, make a blast db:
```
makeblastdb -in db_path/GVDB.fasta -blastdb_version 5 -title "GVDB" -dbtype nucl -out GVDB

```
### 5. Fastq files

You need to specify the path to fastq files of the 10X sequencing.
The script works on 10 Chromium version (2 lanes).

### 6. Config file

Read the config file thoroughly and adjust the variables to your experiment.
The pipeline is saprated to different stages and each stage can run saparatly if the inputs exist.

### 7. Pipeline stages

Each stage can run independently 

#### Database editing

Takes as an input a fasta file of the virus genes you are insteresed to find in your data. It will de-duplicate the database at 99% using CD-HIT and create a cell-ranger reference.

#### Cell ranger count

Uses Cell Ranger to find cells expressing viral genes. Takes as an input fastq files and Cell Ranger reference and return the standard cell-ranger count matrixes.

#### Editing cell ranger count matrix 

Edits the Cell Ranger count matrix: removes empty cells and genes, edits the cell names to include sample names, and saves the output as a pickle file.

#### Choosing infected cells

Chooses infected cells from the Cell Ranger matrix based on three criteria: the total sum of UMIs in a cell, the expression of the most highly-expressed gene,  and the total number of genes expressed in a cell. 
Returns a tab-delimited table of cells and the full path to their related fastq files.

#### Read extraction

Extract the reads of each infected cell from the fastq files and stores them in a cell-specific folder. Reads are trimmed using TrimGalore.
The output is a cell-specific folder with fastq files for each 10X chromium lane.

#### Single-cell assembly

Reads from each cell are assembled using RNA SPAdes. 
The output is a transcript fasta file.

#### Editing reads and assembled transcripts

For each cell, reads from both 10X chromium lanes are combined to a single files and given a cell name (barcode+sample).
Transcript files are edited so the cell name (barcode+sample) will be attached to each transcript name, and all transcript files are combined.
The output is one metatranscriptome fasta file of assembled transcripts from all cells marked with cell names.

#### 18s rRNA Blast

The full metatranscriptome fasta file is compared to the PR2 18s rRNA database.
First, rRNA transcripts are sorted out using sortmerna.
The sorted rRNA transcript file is aligned to the PR2 and metaPR2 databases.
The results are filtered by 99% identity and 100 bp alignment length.
The outputs are Blast results in the 6 format.

#### 18s rRNA lastal

Trimmed reads from each cell are aligned to the PR2 18s rRNA database using lastal.
The outputs are tab-delimited lastal results for each cell.

#### Taxonomic annotation summary

The results from both the Blast and lastal search are combined.
Taxonomic names of specific familes are replaced with the common names (using the replacment file).
Only cells that pass the criteria for both Blast and lastal are considered annotated.
Non-chimeric cells are considered as follows:
1. In the Blast results, the sum of alignments length to one taxa is at least 90% of the total alignment length for a cell.
2. In the lastal results, the sum of reads aligned to one taxa is at least 90% of the total aligned reads for a cell, and the number of aligned reads is at least 10.
3. The results of Blast and lastal are similar.

The outputs are summary of both lastal and blast results, a combined summary of both, and a file of only non-chimeric, annotated cells.

#### Viral Blast

The full metatranscriptome fasta file is compared to a virus database using Blast.
The output is Blast results in the 6 format.

#### Host-virus pairing summary

The results from the viral Blast are summarized and are paired with results from the cell annotation.
The output is a summary of the viral Blast, and a file containing the taxonomic assignment to the virus and to the hosting cell.



