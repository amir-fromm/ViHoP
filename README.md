  
# ViHoP: Virus-Host Pairing
A bioinformatic pipeline for pairing viruses to their hosts from environmental single-cell 10x data.

amir.fromm@weizmann.ac.il  
## What the pipeline does
This pipeline takes single-cell RNA-sequencing data 10X sequencing and, using viral and host references, produces a matrix of host cells and the viruses that infect them.
This repository is built to work on an LSF platform since many of the scripts use parallelization.

![alt text](https://github.com/amir-fromm/ViHoP/blob/main/sc_virus_host_detection_pipeline.png?raw=true)

## 1. Installation
### Cloning this repository

```
git clone https://github.com/amir-fromm/ViHoP.git
``` 
Then, go into the directory.
```
cd ViHoP
```
Create an environment with the requirement.txt file
```
conda create --name host-virus-pairing --file requirements.txt
conda activate host-virus-pairing
```

### Other dependencies  

Besides the packages from the requirements file, the pipeline also uses modules that were installed in our LSF platform and are called using the module load command:
```
CD-HIT/4.6.6
CellRanger/10.0.0
TrimGalore/0.6.5
SPAdes/3.15
BBMap
BLAST/2.11.0
seqtk

```
The specific module names can be changed using the config file. If you don't have these tools installed on your LSF platform, install them on your conda environment.

## 2. Resources
### Sample table fastq file

This file, found in the docs folder, is a tab-delimited file of the different samples used in the analysis and the full path to their correspomnding 10X fastq files:
sample1	/home/user1/fastq/sample1
sample2	/home/user2/fastq/sample2

### Common names replacement file

replacment_file.tsv is a tab-delimited file found in the docs folder, containing four columns:
level,old_name,new_name,color.
The replacement file is used to replace the names of specific taxonomic groups with common names in the cell annotation files.
The first column is the level of the old taxonomic name you want to replace.
The second column is the old taxonomic name.
The third column is the new name you want to assign to the taxonomic group.
The fourth column is a color assigned to the taxonomic group.
For example:

Supergroup	Chlorophyta	Green-algae/Plants	#008000

Result: Cells belonging to the supergroup Chlorophyta will be annotated as Green-algae/Plants and marked green.

### Reference gene set

A reference gene set is required to map the scRNA-seq data to.


### Making host reference databases

You will need to to download the latest protist host 18s rRNA databases PR2 (version 5) and metaPR2 (version 2) (https://pr2-database.org/, https://app.metapr2.org/metapr2/) and create last and blast db from them:
```
makeblastdb -in db_path/PR2.fasta -blastdb_version 5 -title "PR2" -dbtype nucl -out PR2
makeblastdb -in db_path/metaPR2.fasta -blastdb_version 5 -title "metaPR2" -dbtype nucl -out metaPR2
lastdb -c -P 4 PR2 db_path/PR2.fasta

```
### Virus taxonomy

A virus database of your choice should be available (you can use Aylward's giant virus database: https://faylward.github.io/GVDB/), and a virus taxonomy table should be available. 
It is a comma-delimited file where the first column contains the contig name in the virus reference (or part of it), and the second is the virus taxonomy.
After you download the virus database, make a blast db:
```
makeblastdb -in db_path/GVDB.fasta -blastdb_version 5 -title "GVDB" -dbtype nucl -out GVDB

```
### Fastq files

You need to specify the path to fastq files of the 10X sequencing.
The script works on 10x Chromium with 2 lanes.

## 3. Running the pipeline
### Config file

Read the config file thoroughly and adjust the variables to your experiment.
The pipeline is separated into different stages, and each stage can run separately if the inputs exist.

### Pipeline stages

Each stage can run independently if it has the right inputs. You can skip each stage by changing the RUN parameter associated with each stage from "true" to "false". To run the entire pipeline, change all RUN parameters to true. 
In the first run, it is recommended to run the pipeline step by step.

#### Database Editing

Editing the reference gene databse to make it suitable for Cell Ranger count. It takes as an input a fasta file of the virus genes you are interested in finding in your data. It will de-duplicate the database using CD-HIT and create a Cell Ranger reference.

#### Cell ranger count

Using Cell Ranger to find cells expressing viral genes. It takes as an input fastq files and Cell Ranger reference and returns the standard cell-ranger count matrixes.

#### Editing Cell Ranger count matrix 

Editing the Cell Ranger count matrix.  and saves the output as a pickle file.

#### Choosing cells

Concatanate cell ranger matrixes from all samples. Edits the cell names to include sample names.
Cells from the Cell Ranger matrix are chosen based on the total sum of UMIs in a cell. Genes are filtered based on number of cells expressing it.
The output is a tab-delimited table of cells with the full path to their related fastq files.

#### Read extraction

Extracting the reads of each infected cell from the fastq files and storing them in a cell-specific folder. Reads are trimmed using TrimGalore.
The output is a cell-specific folder with fastq files for each 10X chromium lane.

#### Single-cell assembly

Assembling reads from each cell using RNA SPAdes. 
The output is a transcript fasta file.

#### Editing reads and assembled transcripts

For each cell, reads from both 10X chromium lanes are combined into a single file and given a cell name (barcode+sample).
Transcript files are edited so the cell name (barcode+sample) will be attached to each transcript name, and all transcript files are combined.
The output is one metatranscriptome fasta file of assembled transcripts from all cells marked with cell names.

#### 18s rRNA Blast

The combined metatranscriptome fasta file is compared to the PR2 18s rRNA database.
First, rRNA transcripts are sorted out using barrnap.
The sorted rRNA transcript file is aligned to the PR2 and metaPR2 databases.
The outputs are Blast results in the 6 format.

#### 18s rRNA lastal

Trimmed reads from each cell are aligned to the PR2 18s rRNA database using lastal.
The outputs are tab-delimited lastal results for each cell.
The last results are filtered by custom identity and bp alignment length.


#### Taxonomic annotation summary

The Blast results are pre-filtered by custom identity and bp alignment length.
The results from both the Blast and lastal search are combined and summarized according to a taxonomic level.
The best hit for each contig is filtered using a custom threshold and bp alignment length.
If "Annotation" is chosen as a taxonomic level, taxonomic names of specific taxonomic groups are replaced with common names (using the replacement file).
Non-chimeric cells (for pairing) are considered as follows:
1. In the Blast results, no other assembled contig with a different annotation
2. In the lastal results, the sum of reads aligned to one taxa is at least 80% of the total aligned reads for a cell, and the number of aligned reads is at least 10.
3. The results of Blast and lastal are similar.
Only cells that pass the criteria for both Blast and lastal are considered annotated.

The outputs are summaries of lastal and Blast results, a combined summary of both, and a file of only non-chimeric annotated cells.

#### Viral Blast

The combined metatranscriptome fasta file is compared to a virus database using lastal.
The suspected viral contigs are compared to a general reference database (Refseq) using lastal and contigs that matched best to organisms other than viruses are removed ad false positives.
The cleaned tables are merged with best hits retained.
The output is Blast results in the 6 format.

#### Host-virus pairing summary

The results from the viral lastal search are summarized and are paired with results from the cell annotation.
The output is a summary of the viral lastal search with the taxonomic assignment of the virus and the hosting cell.
