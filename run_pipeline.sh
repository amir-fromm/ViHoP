#!/bin/sh
# Must pre-install seqkit
# Documents in docs:
    # docs/sample_table_fastq.tsv # Name of samples and their exact path
    # docs/replacement_file.txt
. ./pairing.config

mkdir -p $WD/log
mkdir -p $WD/tmp

# Generate unique filenames using the current date and time
timestamp=$(date +"%Y%m%d_%H%M%S")
stdout_log="pipeline_log_${timestamp}.o"
stderr_log="pipeline_log_${timestamp}.e"

# Redirect stdout and stderr to the unique log files
exec > >(tee -a "$WD/log/$stdout_log") 2> >(tee -a "$WD/log/$stderr_log" >&2)


if [ "$RUN_DB_EDIT" == "true" ]
then
    echo "Running database editing"
    sleep 1
    eval "$CD_HIT"
    eval "$CELL_RANGER"
    rm -r $REF_NAME
    DB_VIRUS_OUTPUT=$(basename "$DB_VIRUS" .fasta)
    
    if ! ls "$DB_VIRUS_OUTPUT.$DEDUP.fasta" 1> /dev/null 2>&1
    then
        job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/CD-HIT. -o $WD/log/CD-HIT.%J.o -e $WD/log/CD-HIT.%J.e  -n $CDHIT_CORES -R "rusage[mem=$CDHIT_MEM] span[ptile=$CDHIT_CORES]"  cd-hit-est -c 0.$DEDUP -T 0 -M 0 -i $DB_VIRUS -o $DB_VIRUS_OUTPUT.$DEDUP.fasta | awk '{print $2}' | tr -d '<>')
        job_ids+=($job_id)  # Store the job ID
    
        # Wait for all jobs to finish
        for id in "${job_ids[@]}"
        do
            while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
                sleep 10  # Wait for 10 seconds before checking again
            done # Done checking status of jobs 
        done # The jobs are done
    fi

    tr '|=()' '#' < $DB_VIRUS_OUTPUT.$DEDUP.fasta > $DB_VIRUS_OUTPUT.$DEDUP.clean.fasta

    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $DB_VIRUS_OUTPUT.$DEDUP.clean.fasta  >  $DB_VIRUS_OUTPUT.$DEDUP.final.fasta

    sed -i '1{/^[[:space:]]*$/d}' $DB_VIRUS_OUTPUT.$DEDUP.final.fasta


    bin/create_gtf.sh $DB_VIRUS_OUTPUT.$DEDUP.final.fasta

    sed -i 's/\.\t1\t\./.\t+\t./g' 04_$DB_VIRUS_OUTPUT.$DEDUP.final.fasta.gene_intervals.gtf


    bin/create_reference_new.sh $REF_NAME $DB_VIRUS_OUTPUT.$DEDUP.final.fasta 04_$DB_VIRUS_OUTPUT.$DEDUP.final.fasta.gene_intervals.gtf
    echo "Done database editing"
fi


if [ "$RUN_DB_EDIT" == "false" ]
then
    echo "Skipping database editing"
    sleep 1
fi

if [ "$RUN_CELL_RANGER" == "true" ]
then
    echo "Running Cell ranger count protocol"
    sleep 1
    eval "$CELL_RANGER"
    DIRPATH="$(dirname "${REF}")"
    job_ids=()
    echo 
    if [ "$EXPECT_CELLS" != "false" ]
    then
        if [ "$MULTIPLE_FASTQ" == "true" ]
        then
            for dir in $FASTQDIR
            do
                echo $dir
                ID=$(basename $dir)
                echo $ID
                
                rm -r $ID
                # Submit the job and extract the job ID
                job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/cellranger_$ID. -o $WD/log/cellranger_$ID.%J.o -e $WD/log/cellranger_$ID.%J.e  -n $COUNT_CORES -R "rusage[mem=$COUNT_MEM] span[ptile=$COUNT_CORES]" cellranger count --id=$ID --transcriptome=$REF_PATH/$REF_NAME --fastqs=$dir --project=$ID --create-bam=$CREATE_BAM --expect-cells=$EXPECT_CELLS --chemistry=$CHEMISTRY  --nosecondary | awk '{print $2}' | tr -d '<>')
                job_ids+=($job_id)  # Store the job ID
            done
        else
            for dir in $FASTQDIR*/
            do
                echo $dir
                ID=$(basename $dir)
                echo $ID
                
                rm -r $ID
                # Submit the job and extract the job ID
                job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/cellranger_$ID. -o $WD/log/cellranger_$ID.%J.o -e $WD/log/cellranger_$ID.%J.e  -n $COUNT_CORES -R "rusage[mem=$COUNT_MEM] span[ptile=$COUNT_CORES]" cellranger count --id=$ID --transcriptome=$REF_PATH/$REF_NAME --fastqs=$dir/ --project=$ID --create-bam=$CREATE_BAM --expect-cells=$EXPECT_CELLS --chemistry=$CHEMISTRY  --nosecondary | awk '{print $2}' | tr -d '<>')
                job_ids+=($job_id)  # Store the job ID
            done
        fi
        # Wait for all jobs to finish
        for id in "${job_ids[@]}"
        do
            while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
                sleep 10  # Wait for 10 seconds before checking again
            done # Done checking status of jobs 
        done # The jobs are done
    else
        if [ "$MULTIPLE_FASTQ" == "true" ]
        then
            for dir in $FASTQDIR
            do
                echo $dir
                ID=$(basename $dir)
                echo $ID
                
                rm -r $ID
                # Submit the job and extract the job ID
                job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/cellranger_$ID. -o $WD/log/cellranger_$ID.%J.o -e $WD/log/cellranger_$ID.%J.e  -n $COUNT_CORES -R "rusage[mem=$COUNT_MEM] span[ptile=$COUNT_CORES]" cellranger count --id=$ID --transcriptome=$REF_PATH/$REF_NAME --fastqs=$dir --project=$ID --create-bam=$CREATE_BAM --chemistry=$CHEMISTRY  --nosecondary | awk '{print $2}' | tr -d '<>')
                job_ids+=($job_id)  # Store the job ID
            done
        else
            for dir in $FASTQDIR*/
            do
                echo $dir
                ID=$(basename $dir)
                echo $ID
                
                rm -r $ID
                # Submit the job and extract the job ID
                job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/cellranger_$ID. -o $WD/log/cellranger_$ID.%J.o -e $WD/log/cellranger_$ID.%J.e  -n $COUNT_CORES -R "rusage[mem=$COUNT_MEM] span[ptile=$COUNT_CORES]" cellranger count --id=$ID --transcriptome=$REF_PATH/$REF_NAME --fastqs=$dir/ --project=$ID --create-bam=$CREATE_BAM --chemistry=$CHEMISTRY  --nosecondary | awk '{print $2}' | tr -d '<>')
                job_ids+=($job_id)  # Store the job ID
            done
        fi
        # Wait for all jobs to finish
        for id in "${job_ids[@]}"
        do
            while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
                sleep 10  # Wait for 10 seconds before checking again
            done # Done checking status of jobs 
        done # The jobs are done
    fi
    echo "Cell ranger protocol is done"
fi

if [ "$RUN_CELL_RANGER" == "false" ]
then
    echo "Skipping cell ranger protocol"
    sleep 1
fi

# Edit cell ranger files

if [ "$RUN_READ_MTX" == "true" ]
then
    echo "Running cell ranger matrix editing"
    sleep 1
    rm $WD/tmp/read_mtx.txt
    rm $WD/log/raw_UMI_counts_10X_data_*
    job_ids=()
    ITER=1
    for DIR in */outs/filtered_feature_bc_matrix
    do
        if ! ls "$DIR/data_raw.pickle.gz" 1> /dev/null 2>&1
        then
            SAMPLE=${DIR%%/outs/*}
            echo "$SAMPLE";
            # Submit the job and extract the job ID
            job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/raw_UMI_counts_10X_data_$SAMPLE -o $WD/log/raw_UMI_counts_10X_data_$SAMPLE.%J.o -e $WD/log/raw_UMI_counts_10X_data_$SAMPLE.%J.e -n $CORES_MTX -R "rusage[mem=$MEM_MTX] span[ptile=$CORES_MTX]" python bin/00.00.raw_UMI_counts_10X_data.py --data_dir $DIR --file_type mtx | awk '{print $2}' | tr -d '<>')
            job_ids+=($job_id)  # Store the job ID
        fi
    done
    # Wait for all jobs to finish
    for id in "${job_ids[@]}"
    do
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10  # Wait for 10 seconds before checking again
        done # Done checking status of jobs 
    done # The jobs are done
    job_ids=()
    # combine
    grep -l "reaching LSF memory usage limit" $WD/log/raw_UMI_counts_10X_data_* > $WD/tmp/read_mtx.txt

    while [ -s "$WD/tmp/read_mtx.txt" ] && [ "$ITER" -lt 4 ]
    do
        ((ITER++))
        echo "Jobs failed do to memory. Starting iteration $ITER"
        CORES_MTX=$(( ITER*CORES_MTX ))
        for DIR in */outs/filtered_feature_bc_matrix
        do
            if ! ls "$DIR/data_raw.pickle.gz" 1> /dev/null 2>&1
            then
                SAMPLE=${DIR%%/outs/*}
                echo "$SAMPLE";
                # Submit the job and extract the job ID
                job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/raw_UMI_counts_10X_data_$SAMPLE -o $WD/log/raw_UMI_counts_10X_data_$SAMPLE.%J.o -e $WD/log/raw_UMI_counts_10X_data_$SAMPLE.%J.e -n $CORES_MTX -R "rusage[mem=$MEM_MTX] span[ptile=$CORES_MTX]" python bin/00.00.raw_UMI_counts_10X_data.py --data_dir $DIR --file_type mtx | awk '{print $2}' | tr -d '<>')
                job_ids+=($job_id)  # Store the job ID
            fi
        done
         # Wait for all jobs to finish
        for id in "${job_ids[@]}"
        do
            while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
                sleep 10  # Wait for 10 seconds before checking again
            done # Done checking status of jobs 
        done # The jobs are done
        job_ids=()
        grep -l "reaching LSF memory usage limit" $WD/log/raw_UMI_counts_10X_data_* > $WD/tmp/read_mtx.txt
    done
    echo "Done editing mtx files"
    sleep 1


    ##
    echo "Combine pickle files"
    job_id=()
    job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/combine_pickle -o $WD/log/combine_pickle.%J.o -e $WD/log/combine_pickle.%J.e -n $CORES_MTX -R "rusage[mem=$MEM_MTX] span[ptile=$CORES_MTX]" bin/combine_pickle.py --data_dir "$WD/" --list_folders $DIRS_LIST --output_folder Combined_RAW --raw | awk '{print $2}' | tr -d '<>')
    job_ids+=($job_id)  # Store the job ID

    # Wait for all jobs to finish
    for id in "${job_ids[@]}"
    do
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10  # Wait for 10 seconds before checking again
        done # Done checking status of jobs 
    done # The jobs are done
    job_ids=()
    ##
fi

if [ "$RUN_READ_MTX" == "false" ]
then
    echo "Skipping cell ranger matrix editing"
    sleep 1
fi

# Get a list of infected cells
if [ "$RUN_FIND_INFECTED" == "true" ]
then
    echo "Finding infected cells"
    sleep 1

    job_id=()
    job_id=$(bsub -q $QUEUE_SHORT -J $WD/log/choose_cells -o $WD/log/choose_cells.%J.o -e $WD/log/choose_cells.%J.e -n $CORES_MTX -R "rusage[mem=$MEM_MTX] span[ptile=$CORES_MTX]" bin/choose_cells.py --data_dir Combined_RAW/ --sample_table $SAMPLE_TABLE --output_folder Combined_RAW/ --sum $SUM_UMI --count $GENE_COUNT --exp $MAX_UMI | awk '{print $2}' | tr -d '<>')
    job_ids+=($job_id)  # Store the job ID

    # Wait for all jobs to finish
    for id in "${job_ids[@]}"
    do
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10  # Wait for 10 seconds before checking again
        done # Done checking status of jobs 
    done # The jobs are done
    job_ids=()


    echo "Done searching for infected cells"
fi

if [ "$RUN_FIND_INFECTED" == "false" ]
then
    echo "Skipping finding infected cells"
    sleep 1
fi

if [ "$RUN_EXTRACT_READS" == "true" ]
then
    module purge

    echo "Running read extraction protocol"
    sleep 1

    mkdir -p assemble_first_cells/all_cells
    NEW_LIST_CELLS=${LIST_CELLS%.*}_new.tsv
    tr -d '\r' < $LIST_CELLS > $NEW_LIST_CELLS

    NEW_LIST_CELLS=$(realpath "$NEW_LIST_CELLS")

    cd assemble_first_cells/all_cells/
    eval "$SEQTK"
    eval "$TRIM"

    sh $WD/bin/extract_reads_cell_parallel.sh $NEW_LIST_CELLS
    
    cd $WD
    echo "Done extracting reads"
    
fi

if [ "$RUN_EXTRACT_READS" == "false" ]
then
    echo "Skipping read extraction"
    sleep 1
fi

if [ "$RUN_ASSEMBLE_CELLS" == "true" ]
then
    
    echo "Running cell assembly protocol"
    NEW_LIST_CELLS=${LIST_CELLS%.*}_new.tsv
    tr -d '\r' < $LIST_CELLS > $NEW_LIST_CELLS

    NEW_LIST_CELLS=$(realpath "$NEW_LIST_CELLS")

    MAX_JOBS=$(( MAX_MEM / ( MEM_SPADES*CORES_SPADES ) ))
    ITER=1
    sleep 1
    module purge
    eval "$SPADES"
    mkdir -p $WD/log_spades
    rm $WD/log_spades/assemble_cell_*
    $WD/bin/submit_array_spades.sh $NEW_LIST_CELLS $WD/bin/run_spades_job.sh
    cd $WD
    echo "Done assembling cells"
fi

if [ "$RUN_ASSEMBLE_CELLS" == "false" ]
then
    echo "Skipping cell assembly protocol"
    sleep 1
fi

if [ "$RUN_READ_EDITING" == "true" ]
then
    echo "Running read editing protocol"
    sleep 1
    module purge
    # Combine reads to both lanes

    find $WD/assemble_first_cells/all_cells/ -maxdepth 1 -type d ! -name "." | sort > $WD/tmp/folder_list.txt
    mkdir -p $WD/tmp_batches
    split -d -l 2000 $WD/tmp/folder_list.txt $WD/tmp_batches/folder_batch_
    mkdir -p log_read_editing
    job_ids=()
    for batch_file in $WD/tmp_batches/folder_batch_*
    do
        count=$(wc -l < "$batch_file")
        batch_name=$(basename "$batch_file")

        job_id=$(bsub -q $QUEUE_SHORT -J "combine_fastq_${batch_name}[1-$count]" -o "$WD/log_read_editing/log_combine.${batch_name}.%I.o" -e "$WD/log_read_editing/log_combine.${batch_name}.%I.e" -R "rusage[mem=200]" bin/combine_fastqs.sh "$batch_file" | awk '{print $2}' | tr -d '<>')
        job_ids+=($job_id)
    done
    # Wait for each job
    for id in "${job_ids[@]}"; do
        while bjobs -noheader "$id" 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'
        do
            sleep 10
        done
    done
    rm -r $WD/tmp_batches
    rm -r $WD/tmp/folder_list.txt
    if [ "$DELETE_TMP" == "true" ]
    then
        echo "Deleting temporary sequencing files"
        $WD/delete_files.sh $WD/assemble_first_cells/all_cells/
    fi
    echo "Done editing reads"
fi


if [ "$RUN_READ_EDITING" == "false" ]
then
    echo "Skipping read editing"
    sleep 1
fi

if [ "$RUN_TRANSCRIPTS_EDITING" == "true" ]
then
    echo "Running contig editing protocol"
    sleep 1
    module purge
    cd $WD/assemble_first_cells/all_cells/
    # Edit names of transcripts for each cell
    $WD/bin/array_loop_files.sh $WD/assemble_first_cells/all_cells $WD/bin/editing_contigs_array.sh 100

    cat $WD/assemble_first_cells/all_cells/*/assembly/transcripts.edit.fasta > $WD/assemble_first_cells/all_cells.transcripts.edit.fasta
    cd $WD/assemble_first_cells
    # Remove low-complexity transcripts
    eval "$BBMAP"
    bbduk.sh in=all_cells.transcripts.edit.fasta out=all_cells.transcripts.final.fasta outm=low_complexity.fq entropy=$ENTROPY
    if [[ "$DELETE_TMP" == "true" && -s "$WD/assemble_first_cells/all_cells.transcripts.final.fasta" ]]
    then
        echo "Deleting temporary sequencing files"
        rm $WD/assemble_first_cells/all_cells.transcripts.edit.fasta
        rm -r $WD/assemble_first_cells/all_cells/*/assembly/
    fi
    echo "Done editing contigs"
    
fi

if [ "$RUN_TRANSCRIPTS_EDITING" == "false" ]
then
    echo "Skipping contig editing"
    sleep 1
fi

if [ "$RUN_BLAST" == "true" ]
then
    cd $WD
    echo "Running blast"
    sleep 1
    rm $WD/log/blastn_*pr2*
    eval "$BLAST"
    if [ "$RUN_SORTMERNA" == "true" ]
    then
        mkdir -p assemble_first_cells/sortmerna
        rm -r $WD/kvdb
        sortmerna --ref  $PR2_FASTA --kvdb $WD/kvdb --threads $THREADS_SORTMERNA --reads assemble_first_cells/all_cells.transcripts.final.fasta --aligned assemble_first_cells/sortmerna/all_cells.rRNA -fastx assemble_first_cells/sortmerna/all_cells.rRNA
    fi

    TOTAL_MEM=$(( CORES_BLAST*MEM_BLAST ))
    if [ "$TOTAL_MEM" -lt 2000 ]
    then
        echo "Blast memory has to be higher than 2000!"
        exit 1
    fi
    mkdir -p blast_results
    ITER=1
    THREADS=$(( CORES_BLAST*2 ))
    echo "Starting blast search"
    job_ids=()
    job_id1=$( bsub -q $QUEUE_SHORT -J $WD/log/blastn_pr2   -o $WD/log/blastn_pr2.%J.o -e $WD/log/blastn_pr2.%J.e -n $CORES_BLAST -R "rusage[mem=$MEM_BLAST] span[ptile=$CORES_BLAST]" blastn -outfmt 6 -evalue $EVALUE -query assemble_first_cells/sortmerna/all_cells.rRNA.fa -db  $PR2_DB -num_threads $THREADS -out blast_results/all_cells.rRNA.PR2.tsv | awk '{print $2}' | tr -d '<>')
    job_id2=$( bsub -q $QUEUE_SHORT -J $WD/log/blastn_metapr2 -o $WD/log/blastn_metapr2.%J.o -e $WD/log/blastn_metapr2.%J.e -n $CORES_BLAST -R "rusage[mem=$MEM_BLAST] span[ptile=$CORES_BLAST]" blastn -outfmt 6 -evalue $EVALUE -query assemble_first_cells/sortmerna/all_cells.rRNA.fa -db  $METAPR2_DB -num_threads $THREADS -out blast_results/all_cells.rRNA.metaPR2.tsv | awk '{print $2}' | tr -d '<>')

    job_ids+=($job_id1)  # Store the job ID
    job_ids+=($job_id2)  # Store the job ID

    # Wait for all jobs to finish
    for id in "${job_ids[@]}"; do
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10  # Wait for 10 seconds before checking again
        done # Done checking status of jobs 
    done # The jobs are done
    
    # Check if all some jobs did not have enough memory. If so, increase memory in a loop until all jobs had enough memory.
    grep -l "reaching LSF memory usage limit" $WD/log/blastn_*pr2* > $WD/tmp/blastn_limit_use.txt

    while [ -s "$WD/tmp/blastn_limit_use.txt" ] && [ "$ITER" -lt 6 ]
    do
        ((ITER++))
        MEM_BLAST=$(( MEM_BLAST*ITER ))
        echo "Some jobs failed due to memory. Starting iteration $ITER"
        job_ids=()
            
        job_id1=$( bsub -q $QUEUE_SHORT -J $WD/log/blastn_pr2.$ITER.   -o $WD/log/blastn_pr2.$ITER.%J.o -e $WD/log/blastn_pr2.$ITER.%J.e -n $CORES_BLAST -R "rusage[mem=$MEM_BLAST] span[ptile=$CORES_BLAST]" blastn -outfmt 6 -evalue $EVALUE -query assemble_first_cells/sortmerna/all_cells.rRNA.fa -db  $PR2_DB -num_threads $THREADS -out blast_results/all_cells.rRNA.PR2.tsv | awk '{print $2}' | tr -d '<>')
        job_id2=$( bsub -q $QUEUE_SHORT -J $WD/log/blastn_metapr2.$ITER. -o $WD/log/blastn_metapr2.$ITER.%J.o -e $WD/log/blastn_metapr2.$ITER.%J.e -n $CORES_BLAST -R "rusage[mem=$MEM_BLAST] span[ptile=$CORES_BLAST]" blastn -outfmt 6 -evalue $EVALUE -query assemble_first_cells/sortmerna/all_cells.rRNA.fa -db  $METAPR2_DB -num_threads $THREADS -out blast_results/all_cells.rRNA.metaPR2.tsv | awk '{print $2}' | tr -d '<>')

        job_ids+=($job_id1)  # Store the job ID
        job_ids+=($job_id2)  # Store the job ID

        # Wait for all jobs to finish
        for id in "${job_ids[@]}"
        do
            while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
                sleep 10  # Wait for 10 seconds before checking again
            done # Done checking status of jobs 
        done # The jobs are done
        job_ids=()
        grep -l "reaching LSF memory usage limit" $WD/log/blastn_limit_use*.$ITER.*.o > $WD/tmp/blastn_limit_use.txt
    done
    
    awk '$4 >= 100 && $3 >= 99' blast_results/all_cells.rRNA.PR2.tsv > blast_results/filtered.rRNA.PR2.tsv

    awk '$4 >= 100 && $3 >= 99' blast_results/all_cells.rRNA.metaPR2.tsv > blast_results/filtered.rRNA.metaPR2.tsv
    echo "Done blast"
fi

if [ "$RUN_BLAST" == "false" ]
then
    echo "Skipping blast"
    sleep 1

fi

if [ "$RUN_18S_LAST" == "true" ]
then
    cd $WD
    echo "Running 18s last search"
    sleep 1
    mkdir -p $WD/mapped_cells_18s
    mkdir -p $WD/log_last
    rm $WD/log_last/*
    eval "$SEQTK"
    ls assemble_first_cells/all_cells/*/*combined.fq.gz > $WD/tmp/combined_reads.txt
  #  sh bin/lastal_loop_18s.sh $WD/tmp/combined_reads.txt
    sh bin/lastal_array_18s.sh $WD/tmp/combined_reads.txt
    echo "Done running 18s last search"
fi


if [ "$RUN_18S_LAST" == "false" ]
then
    echo "Skipping 18s last"
    sleep 1
fi


if [ "$RUN_CELL_ANNOTATION" == "true" ]
then
    echo "Running cell annotation"
    sleep 1

    BLAST_OUTPUT=$WD/blast_results/all_cells.rRNA.PR2.tsv

    SZ=$(stat -c '%s' $BLAST_OUTPUT)

    if ((SZ < 100000000))
    then
        echo "file is smaller than 100MB"
        MEM_SUMMARY_BLAST=1000
        MEM_SUMMARY_LAST=4000
        CORES_SUMMARY=5
    elif ((SZ >= 100000000 && SZ <= 1000000000))
    then
        echo "file is between 100MB and 1GB"
        MEM_SUMMARY_BLAST=2000
        MEM_SUMMARY_LAST=8000
        CORES_SUMMARY=5

    elif ((SZ >= 1000000000 && SZ <= 10000000000))
    then
        echo "file is between 1GB and 10GB"
        MEM_SUMMARY_BLAST=4000
        MEM_SUMARY_LAST=16000
        CORES_SUMMARY=5

    else
        echo "file is larger than 10GB"
        MEM_SUMMARY_BLAST=5000
        MEM_SUMARY_LAST=13000
        CORES_SUMMARY=12
    fi
    job_ids=()
    if [ "$SUMMARIZE_BLAST" == "true" ]
    then
        $WD/bin/tsv_to_parquet.py --tsv_path $WD/blast_results/filtered.rRNA.metaPR2.tsv
        $WD/bin/tsv_to_parquet.py --tsv_path $WD/blast_results/filtered.rRNA.PR2.tsv

        job_id1=$( bsub -q $QUEUE_SHORT -J $WD/log/blast_summary.   -o $WD/log/blast_summary.%J.o -e $WD/log/blast_summary.%J.e -n $CORES_SUMMARY -R "rusage[mem=$MEM_SUMMARY_BLAST] span[ptile=$CORES_SUMMARY]" bin/blast_sum_parquet.py --path $WD --metapr2_file /blast_results/filtered.rRNA.metaPR2.parquet --pr2_file /blast_results/filtered.rRNA.PR2.parquet --replacement_file $REPLACEMENT_FILE --identity $BLAST_IDENTITY --length $BLAST_LENGTH --level $TAXONOMIC_LEVEL | awk '{print $2}' | tr -d '<>')
        job_ids+=($job_id1)  # Store the job ID
    fi
    if [ "$SUMMARIZE_LAST" == "true" ]
    then
        job_id2=$( bsub -q $QUEUE_SHORT -J $WD/log/last_summary. -o $WD/log/last_summary.%J.o -e $WD/log/last_summary.%J.e -n $CORES_SUMMARY -R "rusage[mem=$MEM_SUMARY_LAST] span[ptile=$CORES_SUMMARY]" bin/last_sum_18s_multi.py --last_path "mapped_cells_18s/" --last_summary $LAST_SUMMARY --replacement_file $REPLACEMENT_FILE --identity $LAST_IDENTITY --length $LAST_LENGTH --level $TAXONOMIC_LEVEL --workers $CORES_SUMMARY | awk '{print $2}' | tr -d '<>')
        job_ids+=($job_id2)  # Store the job ID
    fi
    # Wait for all jobs to finish
    for id in "${job_ids[@]}"; do
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10  # Wait for 10 seconds before checking again
        done # Done checking status of jobs 
    done # The jobs are done
    PATH_LAST_SUM=last_summary."$TAXONOMIC_LEVEL".csv 
    PATH_BLAST_SUM=blast_summary."$TAXONOMIC_LEVEL".csv 
    LAST_COUNT_FINAL=$(( LAST_COUNT*2 )) # The minimum number of reads mapped in last must be at least double the lower threshold (to avoid the case of 10 reads mapped to one, 9 mapped to second )
    bin/summarize_annotations.py --last_path $PATH_LAST_SUM --blast_path $PATH_BLAST_SUM --cell_path $LIST_CELLS --count_reads $LAST_COUNT_FINAL --percent $THRESHOLD_CHIMERIC --level $TAXONOMIC_LEVEL
    if [ -e "$PATH_LAST_SUM" ] && [ "$DELETE_TMP" == "true" ]
    then
        rm 18s_last_df.csv
        rm -r mapped_cells_18s/
    fi 
    echo "Done annotating cells"
fi

if [ "$RUN_CELL_ANNOTATION" == "false" ]
then
    echo "Skipping cell annotation"
    sleep 1
fi

if [ "$RUN_VIRAL_BLAST" == "true" ]
then
    cd $WD
    echo "Running viral blast"
    rm $WD/log/blastn_virus*
    sleep 1
    eval "$BLAST"

    if ! ls "$BLAST_DB_VIRUS.nhr" 1> /dev/null 2>&1 ; then # check if s specific database file exists
        DIR="$(dirname "${BLAST_DB_VIRUS}")"
        FILE="$(basename "${BLAST_DB_VIRUS}")"
        echo "Creating blast database for $FILE in directory $DIR"
        makeblastdb -in $BLAST_DB_VIRUS_FASTA -blastdb_version 5 -title "$FILE" -dbtype nucl -out $BLAST_DB_VIRUS
    fi
    mkdir -p blast_results
    ITER=1
    THREADS=$(( CORES_BLAST*2 ))
    
    job_ids=()
    job_id=$( bsub -q $QUEUE_SHORT -J $WD/log/blastn_virus.   -o $WD/log/blastn_virus.%J.o -e $WD/log/blastn_virus.%J.e -n $CORES_BLAST -R "rusage[mem=$MEM_BLAST] span[ptile=$CORES_BLAST]" blastn -outfmt 6 -word_size 11 -evalue $EVALUE -query assemble_first_cells/all_cells.transcripts.final.fasta -db  $BLAST_DB_VIRUS -num_threads $THREADS -out blast_results/all_cells.virus.tsv | awk '{print $2}' | tr -d '<>')

    job_ids+=($job_id)  # Store the job ID

    # Wait for all jobs to finish
    for id in "${job_ids[@]}"; do
        while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
            sleep 10  # Wait for 10 seconds before checking again
        done # Done checking status of jobs 
    done # The jobs are done
    
    # Check if all some jobs did not have enough memory. If so, increase memory in a loop until all jobs had enough memory.
    
    grep -l "reaching LSF memory usage limit" $WD/log/blastn_virus* > $WD/tmp/blastn_limit_use.txt

    while [ -s "$WD/tmp/blastn_limit_use.txt" ] && [ "$ITER" -lt 6 ]
    do
        ((ITER++))
        MEM_BLAST=$(( MEM_BLAST*ITER ))
        echo "Some jobs failed due to memory. Starting iteration $ITER"
        job_ids=()
            
        job_id=$( bsub -q $QUEUE_SHORT -J $WD/log/blastn_virus.   -o $WD/log/blastn_virus.%J.o -e $WD/log/blastn_virus.%J.e -n $CORES_BLAST -R "rusage[mem=$MEM_BLAST] span[ptile=$CORES_BLAST]" blastn -outfmt 6 -evalue $EVALUE -word_size 11 -query assemble_first_cells/all_cells.transcripts.final.fasta -db  $BLAST_DB_VIRUS -num_threads $THREADS -out blast_results/all_cells.virus.tsv | awk '{print $2}' | tr -d '<>')
        job_ids+=($job_id)  # Store the job ID

        # Wait for all jobs to finish
        for id in "${job_ids[@]}"
        do
            while bjobs -noheader $id 2>/dev/null | grep -q 'PEND\|RUN\|UNKWN\|SSUSP'; do
                sleep 10  # Wait for 10 seconds before checking again
            done # Done checking status of jobs 
        done # The jobs are done
        job_ids=()
        grep -l "reaching LSF memory usage limit" $WD/log/blastn_limit_use*.$ITER.*.o > $WD/tmp/blastn_limit_use.txt
    done
    echo "Done virus blast"

fi
    

if [ "$RUN_VIRAL_BLAST" == "false" ]
then
    echo "Skipping viral blast"
    sleep 1
fi


if [ "$RUN_VIRAL_SUMMARY" == "true" ]
then
    sleep 1
    echo "Running viral summary"
    bin/blast_viruses_summary.py --output_path $WD --output_name summary_blast_viruses.$TAXONOMIC_LEVEL.tsv --taxa_file $VIRUS_TAXONOMY --blast_file $VIRUS_BLAST_FILE --annotation_file summarized_annotation.$TAXONOMIC_LEVEL.csv --identity $VIRUS_IDENTITY --length $VIRUS_LENGTH
    echo "Done viral summary"
fi

if [ "$RUN_VIRAL_SUMMARY" == "false" ]
then
    echo "Skipping viral summary"
    sleep 1
fi
